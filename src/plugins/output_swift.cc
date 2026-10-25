// This file is part of MUSIC2, a code for generating multi-scale
// cosmological initial conditions.
//
// Copyright (C) 2010 Oliver Hahn
// Copyright (C) 2026 Will Roper (SWIFT cell metadata and provenance)
//
// MUSIC2 is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.

#ifdef HAVE_HDF5

#define GAS_PARTTYPE 0
#define HIGHRES_DM_PARTTYPE 1
#define COARSE_DM_DEFAULT_PARTTYPE 2
#define STAR_PARTTYPE 4
#define NTYPES 7

#include "HDF_IO.hh"
#include "cmake_config.hh"
#include "cosmology_calculator.hh"
#include "output.hh"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <sstream>
#include <string>
#ifdef USE_MPI
#include <mpi.h>
#endif

// The legacy output-plugin factory passes config only; main owns initialized
// cosmology calculator used here to persist the effective (defaulted) values.
extern std::unique_ptr<cosmology::calculator> the_cosmo_calc;

/**
 * @brief Write SWIFT-compatible HDF5 initial conditions.
 *
 * MUSIC2 streams grid-derived particle fields one component at a time. This
 * plugin buffers each field, partitions the replicated arrays across MPI ranks
 * at finalize(), then applies a stable per-rank cell permutation. Global
 * per-cell counts and rank offsets produce cell-major, rank-major datasets and
 * the /Cells lookup table used by SWIFT and swiftsimio.
 *
 * Parallel HDF5 writes collectively to the shared file. With serial HDF5,
 * ranks take turns writing disjoint cell hyperslabs. Both modes produce one
 * file and use whole-box geometry for every particle type.
 */
class swift_output_plugin : public output_plugin {
protected:
  // header/config
  std::vector<unsigned long long> nPart;
  std::vector<long long> nPartTotal;
  std::vector<double> massTable;
  double time, redshift, boxSize;
  unsigned int numFiles, coarsePartType;

  double omega0, omega_L, hubbleParam;

  // configuration
  double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
  double omega_b, rhoCrit;
  double posFac, velFac;
  long long nPartTotAllTypes;
  bool doBaryons, useLongIDs;

  //! Number of SWIFT top-level cells along each dimension.
  size_t cdim_;

  //! Total number of top-level cells (cdim_^3).
  size_t ncells_;

  //! Box side length in output length units.
  double dim_;

  //! Top-level cell side length in output length units.
  double cell_width_;

  //! Inverse top-level cell width, matching SWIFT's iwidth.
  double iwidth_;

  //! Full replicated particle fields collected from MUSIC's component-wise API.
  struct particle_buffer {
    //! Flat (particle, xyz) arrays and scalar fields in MUSIC's original order.
    std::vector<double> coordinates, velocities, masses, smoothing_lengths,
        internal_energies;
    //! IDs are held wide until the selected on-disk type is written.
    std::vector<unsigned long long> ids;
  } buffers_[NTYPES];

  //! MPI rank layout; {0,1} in a non-MPI build.
  int rank_, nranks_;

  size_t npfine, npart, npcoarse;
  size_t npfine_dm, npfine_gas; // finest-level particles per species (differ from npfine for masked loads)
  std::vector<size_t> levelcounts;

  //! Masked particle load ('ParticleLoad = masked'); disabled means one DM and one gas particle per fine cell.
  particle_load_mask pmask;

  // parameter file hints
  int pmgrid, gridboost;
  float softening, Tini;
  double gamma, YHe, Tcmb0; //for gas properties

  using output_plugin::cf_;

  /**
   * @brief Attach SWIFT unit exponents to one particle dataset.
   *
   * Current, temperature, and Hubble-factor exponents are zero for every
   * field written by this plugin.
   */
  void write_unit_attributes(const std::string &filename, const std::string &dataset,
                             double mass, double length, double timeExponent,
                             double scaleFactor) {
    HDFWriteDatasetAttribute(filename, dataset, "U_M exponent", std::vector<double>(1, mass));
    HDFWriteDatasetAttribute(filename, dataset, "U_L exponent", std::vector<double>(1, length));
    HDFWriteDatasetAttribute(filename, dataset, "U_t exponent", std::vector<double>(1, timeExponent));
    HDFWriteDatasetAttribute(filename, dataset, "U_I exponent", std::vector<double>(1, 0.0));
    HDFWriteDatasetAttribute(filename, dataset, "U_T exponent", std::vector<double>(1, 0.0));
    HDFWriteDatasetAttribute(filename, dataset, "a-scale exponent", std::vector<double>(1, scaleFactor));
    HDFWriteDatasetAttribute(filename, dataset, "h-scale exponent", std::vector<double>(1, 0.0));
  }

  /**
   * @brief Write the whole-box SWIFT top-level cell tree.
   *
   * Centres and bounds are geometric cell positions, not particle extrema.
   * This lets swiftsimio select cells without padding a requested region.
   * Files is zero everywhere because this implementation writes one file.
   */
  void create_cell_tree(const std::string &filename) {
    std::vector<double> centres(3 * ncells_), minima(3 * ncells_), maxima(3 * ncells_);
    for (size_t i = 0; i < cdim_; ++i) {
      for (size_t j = 0; j < cdim_; ++j) {
        for (size_t k = 0; k < cdim_; ++k) {
          const size_t cell = k + cdim_ * (j + cdim_ * i);
          const size_t index[3] = {i, j, k};
          for (size_t d = 0; d < 3; ++d) {
            minima[3 * cell + d] = index[d] * cell_width_;
            centres[3 * cell + d] = (index[d] + 0.5) * cell_width_;
            maxima[3 * cell + d] = (index[d] + 1.0) * cell_width_;
          }
        }
      }
    }

    HDFCreateGroup(filename, "Cells");
    HDFCreateSubGroup(filename, "Cells", "Meta-data");
    HDFCreateSubGroup(filename, "Cells", "Files");
    HDFCreateSubGroup(filename, "Cells", "OffsetsInFile");
    HDFCreateSubGroup(filename, "Cells", "Counts");
    HDFCreateSubGroup(filename, "Cells", "MinPositions");
    HDFCreateSubGroup(filename, "Cells", "MaxPositions");
    HDFWriteGroupAttribute(filename, "/Cells/Meta-data", "nr_cells",
                           std::vector<int>(1, static_cast<int>(ncells_)));
    HDFWriteGroupAttribute(filename, "/Cells/Meta-data", "size", std::vector<double>(3, cell_width_));
    HDFWriteGroupAttribute(filename, "/Cells/Meta-data", "dimension", std::vector<int>(3, cdim_));
    HDFWriteGroupAttribute(filename, "/Cells/Meta-data", "Origin", std::vector<double>(3, 0.0));
    HDFWriteDatasetVector(filename, "/Cells/Centres", centres);

    for (unsigned type = 0; type < NTYPES; ++type) {
      const particle_buffer &buffer = buffers_[type];
      if (buffer.coordinates.empty()) continue;
      const std::string typeName = "PartType" + std::to_string(type);
      HDFWriteDatasetVector(filename, "/Cells/MinPositions/" + typeName, minima);
      HDFWriteDatasetVector(filename, "/Cells/MaxPositions/" + typeName, maxima);
    }
  }

  /**
   * @brief Store reproducibility metadata under /ICs_parameters.
   *
   * SWIFT preserves this group in later snapshots. Config contains every
   * stored option, one subgroup per section; Cosmology records effective
   * parameter values; Build records compiler and library versions; Derived
   * records particle/gas values; input spectra and CLASS parameters are copied
   * when MUSIC generated corresponding sidecar files.
   */
  void write_provenance_metadata(const std::string &filename) {
    // Preserve every stored configuration option exactly as MUSIC sees it.
    HDFCreateSubGroup(filename, "ICs_parameters", "Config");
    std::string previousSection;
    for (const auto &item : cf_.get_all()) {
      const size_t separator = item.first.find('/');
      const std::string section = item.first.substr(0, separator);
      const std::string key = item.first.substr(separator + 1);
      if (section.empty() || key.empty()) {
        continue;
      }
      if (section != previousSection) {
        HDFCreateSubGroup(filename, "/ICs_parameters/Config", section);
        previousSection = section;
      }
      HDFWriteGroupAttribute(filename, "/ICs_parameters/Config/" + section, key, item.second);
    }

    // Record the effective cosmological values consumed by this plugin.
    HDFCreateSubGroup(filename, "ICs_parameters", "Cosmology");
    for (const auto &parameter : the_cosmo_calc->cosmo_param_.get_all()) {
      HDFWriteGroupAttribute(filename, "/ICs_parameters/Cosmology",
                             parameter.first, parameter.second);
    }

    // Record enough build information to identify numerical/library changes.
    HDFCreateSubGroup(filename, "ICs_parameters", "Build");
    HDFWriteGroupAttribute(filename, "/ICs_parameters/Build", "Build Type", std::string(CMAKE_BUILDTYPE_STR));
    HDFWriteGroupAttribute(filename, "/ICs_parameters/Build", "Compiler", std::string(__VERSION__));
    HDFWriteGroupAttribute(filename, "/ICs_parameters/Build", "FFTW Version",
                           std::string(FFTW_API(version)));
    unsigned hdf5Major, hdf5Minor, hdf5Release;
    H5get_libversion(&hdf5Major, &hdf5Minor, &hdf5Release);
    HDFWriteGroupAttribute(filename, "/ICs_parameters/Build", "HDF5 Version",
                           std::to_string(hdf5Major) + "." + std::to_string(hdf5Minor) + "." +
                               std::to_string(hdf5Release));
    HDFWriteGroupAttribute(filename, "/ICs_parameters/Build", "MPI Ranks", nranks_);
    HDFWriteGroupAttribute(filename, "/ICs_parameters/Build", "Threads",
                           CONFIG::num_threads);
#ifdef USE_PARALLEL_HDF5
    HDFWriteGroupAttribute(filename, "/ICs_parameters/Build", "HDF5 I/O",
                           std::string("parallel collective"));
#else
    HDFWriteGroupAttribute(filename, "/ICs_parameters/Build", "HDF5 I/O",
                           std::string("serialized serial"));
#endif

    HDFCreateSubGroup(filename, "ICs_parameters", "Derived");
    const double targetRedshift =
        cf_.get_value_safe<double>("cosmology", "ztarget", 0.0);
    HDFWriteGroupAttribute(filename, "/ICs_parameters/Derived", "D+ start",
                           the_cosmo_calc->get_growth_factor(time));
    HDFWriteGroupAttribute(
        filename, "/ICs_parameters/Derived", "D+ target",
        the_cosmo_calc->get_growth_factor(1.0 / (1.0 + targetRedshift)));
    HDFWriteGroupAttribute(
        filename, "/ICs_parameters/Derived", "sigma_8",
        static_cast<double>(the_cosmo_calc->compute_sigma8()));
    for (unsigned type = 0; type < NTYPES; ++type) {
      const auto &masses = buffers_[type].masses;
      if (masses.empty()) continue;
      const double mean = std::accumulate(masses.begin(), masses.end(), 0.0) / masses.size();
      HDFWriteGroupAttribute(filename, "/ICs_parameters/Derived",
                             "Mean particle mass PartType" + std::to_string(type) + " [U_M]", mean);
    }
    if (doBaryons) {
      const double mu = (Tini > 1.e4) ? 4.0 / (8. - 5. * YHe) : 4.0 / (1. + 3. * (1. - YHe));
      const double npol = (std::fabs(1.0 - gamma) > 1e-7) ? 1.0 / (gamma - 1.0) : 1.0;
      const double energy = 1.380649e-16 / 1.67262192369e-24 * Tini * npol / mu /
                            UnitVelocity_in_cm_per_s / UnitVelocity_in_cm_per_s;
      HDFWriteGroupAttribute(filename, "/ICs_parameters/Derived", "Gas temperature over mu [K]", Tini / mu);
      HDFWriteGroupAttribute(filename, "/ICs_parameters/Derived",
                             "Gas internal energy [(U_L/U_t)^2]", energy);
      HDFWriteGroupAttribute(filename, "/ICs_parameters/Derived", "Gas smoothing length [U_L]",
                             gas_smoothing_length());
    }

    std::ifstream powerspectrum("input_powerspec.txt");
    if (powerspectrum) {
      std::vector<std::string> names;
      std::vector<std::vector<double>> columns;
      std::string line;
      while (std::getline(powerspectrum, line)) {
        if (!line.empty() && line[0] == '#') {
          if (line.find("k [h/Mpc]") != std::string::npos) {
            std::istringstream header(line.substr(1));
            std::string name;
            while (header >> name) {
              if (name == "[h/Mpc]") continue;
              names.push_back(name);
            }
          }
          continue;
        }
        std::istringstream values(line);
        double value;
        size_t column = 0;
        while (values >> value) {
          if (columns.size() <= column) columns.emplace_back();
          columns[column++].push_back(value);
        }
      }
      if (!columns.empty()) {
        HDFCreateSubGroup(filename, "ICs_parameters", "InputPowerSpectrum");
        HDFWriteGroupAttribute(filename, "/ICs_parameters/InputPowerSpectrum", "k units",
                               std::string("h/Mpc"));
        HDFWriteGroupAttribute(filename, "/ICs_parameters/InputPowerSpectrum", "Scale-factor", time);
        for (size_t c = 0; c < columns.size(); ++c) {
          const std::string name = c < names.size() ? names[c] : "column_" + std::to_string(c);
          HDFWriteDataset(filename, "/ICs_parameters/InputPowerSpectrum/" + name, columns[c]);
        }
      }
    }

    // The CLASS transfer-function plugin writes its complete input key/value
    // list alongside the ICs. Preserve it verbatim as string attributes.
    if (cf_.get_value_safe<std::string>("cosmology", "transfer", "") ==
        "CLASS") {
      std::ifstream classInput("input_class_parameters.ini");
      if (classInput) {
        HDFCreateSubGroup(filename, "ICs_parameters", "CLASS");
        std::string line;
        while (std::getline(classInput, line)) {
          const size_t equals = line.find('=');
          if (equals == std::string::npos)
            continue;
          auto trim = [](std::string value) {
            const size_t first = value.find_first_not_of(" \t\r");
            const size_t last = value.find_last_not_of(" \t\r");
            return first == std::string::npos
                       ? std::string()
                       : value.substr(first, last - first + 1);
          };
          const std::string key = trim(line.substr(0, equals));
          const std::string value = trim(line.substr(equals + 1));
          if (!key.empty() && !value.empty())
            HDFWriteGroupAttribute(filename, "/ICs_parameters/CLASS", key,
                                   value);
        }
      }
      HDFWriteGroupAttribute(
          filename, "ICs_parameters", "Target Redshift",
          cf_.get_value_safe<double>("cosmology", "ztarget", 0.0));
    }

    // PANPHASIA realizations are identified by their descriptor.
    if (cf_.get_value_safe<std::string>("random", "generator", "") ==
        "PANPHASIA") {
      HDFWriteGroupAttribute(
          filename, "ICs_parameters", "Descriptor",
          cf_.get_value<std::string>("random", "descriptor"));
    }
  }

  /**
   * @brief Create one full-size particle dataset before rank writes begin.
   *
   * Serial HDF5 uses bounded chunks and Fletcher32 checksums. Parallel HDF5
   * keeps datasets unfiltered so every rank can participate in collective I/O.
   */
  template <typename T>
  void create_particle_dataset(const std::string &path, size_t particles,
                               size_t width = 1) {
    hid_t file = H5Fopen(fname_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    hid_t type = GetDataType<T>();
    hid_t properties = H5Pcreate(H5P_DATASET_CREATE);
    hid_t space;
    if (width == 1) {
      hsize_t dims[1] = {particles};
      space = H5Screate_simple(1, dims, NULL);
#ifndef USE_PARALLEL_HDF5
      hsize_t chunk[1] = {std::max<hsize_t>(1, std::min<hsize_t>(1024 * 1024 / sizeof(T), particles))};
      H5Pset_chunk(properties, 1, chunk);
      H5Pset_fletcher32(properties);
#endif
    } else {
      hsize_t dims[2] = {particles, width};
      space = H5Screate_simple(2, dims, NULL);
#ifndef USE_PARALLEL_HDF5
      hsize_t chunk[2] = {std::max<hsize_t>(1, std::min<hsize_t>(256 * 1024 / sizeof(T), particles)), width};
      H5Pset_chunk(properties, 2, chunk);
      H5Pset_fletcher32(properties);
#endif
    }
    hid_t dataset = H5Dcreate(file, path.c_str(), type, space, properties);
    if (dataset < 0) {
      H5Sclose(space);
      H5Pclose(properties);
      H5Fclose(file);
      throw std::runtime_error("Failed to create SWIFT dataset " + path);
    }
    H5Dclose(dataset);
    H5Sclose(space);
    H5Pclose(properties);
    H5Fclose(file);
  }

  /**
   * @brief Write one field into this rank's disjoint cell hyperslabs.
   *
   * localCounts and localOffsets describe the rank's packed cell-sorted buffer;
   * rankOffsets give the corresponding global file starts. OR selections let
   * one H5Dwrite cover all populated cells, including empty selections for
   * ranks with no particles in this species.
   */
  template <typename T>
  void write_cell_ordered_dataset(hid_t file, const std::string &path,
                                  const std::vector<T> &values,
                                  const std::vector<long long> &localCounts,
                                  const std::vector<long long> &rankOffsets,
                                  const std::vector<long long> &localOffsets,
                                  size_t width, hid_t dxpl) {
    hid_t dataset = H5Dopen(file, path.c_str());
    hid_t fileSpace = H5Dget_space(dataset);
    hid_t memorySpace;
    if (values.empty()) {
      memorySpace = H5Screate(H5S_NULL);
      H5Sselect_none(fileSpace);
    } else {
      hsize_t memoryDims[2] = {values.size() / width, width};
      memorySpace = H5Screate_simple(width == 1 ? 1 : 2, memoryDims, NULL);
      bool first = true;
      for (size_t cell = 0; cell < ncells_; ++cell) {
        if (!localCounts[cell]) continue;
        hsize_t fileStart[2] = {static_cast<hsize_t>(rankOffsets[cell]), 0};
        hsize_t memoryStart[2] = {static_cast<hsize_t>(localOffsets[cell]), 0};
        hsize_t count[2] = {static_cast<hsize_t>(localCounts[cell]), width};
        H5Sselect_hyperslab(fileSpace, first ? H5S_SELECT_SET : H5S_SELECT_OR,
                            fileStart, NULL, count, NULL);
        H5Sselect_hyperslab(memorySpace, first ? H5S_SELECT_SET : H5S_SELECT_OR,
                            memoryStart, NULL, count, NULL);
        first = false;
      }
      if (first) {
        H5Sselect_none(fileSpace);
        H5Sselect_none(memorySpace);
      }
    }
    if (H5Dwrite(dataset, GetDataType<T>(), memorySpace, fileSpace, dxpl,
                 values.empty() ? NULL : values.data()) < 0) {
      H5Sclose(memorySpace);
      H5Sclose(fileSpace);
      H5Dclose(dataset);
      throw std::runtime_error("Failed to write SWIFT dataset " + path);
    }
    H5Sclose(memorySpace);
    H5Sclose(fileSpace);
    H5Dclose(dataset);
  }

  /**
   * @brief Partition, stable-sort, and write one SWIFT particle type.
   *
   * Each rank receives a contiguous slice of MUSIC's replicated global
   * particle list. Per-cell Allreduce and Exscan results establish global
   * cell offsets and rank-specific file ranges. Sorting is stable inside each
   * rank, so resulting file order is cell-major then rank-major.
   */
  void write_particle_type(unsigned type) {
    particle_buffer &buffer = buffers_[type];
    const size_t total = buffer.coordinates.size() / 3;
    if (!total) return;
    if (buffer.coordinates.size() % 3 != 0 ||
        buffer.velocities.size() != 3 * total ||
        buffer.masses.size() != total || buffer.ids.size() != total ||
        (type == GAS_PARTTYPE && doBaryons &&
         (buffer.smoothing_lengths.size() != total ||
          buffer.internal_energies.size() != total))) {
      throw std::runtime_error("Incomplete SWIFT particle fields for PartType" +
                               std::to_string(type));
    }
    const size_t begin = total * static_cast<size_t>(rank_) / nranks_;
    const size_t end = total * static_cast<size_t>(rank_ + 1) / nranks_;
    const size_t local = end - begin;
    std::vector<size_t> cellIds(local), permutation(local);
    std::vector<long long> localCounts(ncells_, 0), globalCounts(ncells_, 0),
        preceding(ncells_, 0), cellOffsets(ncells_, 0), rankOffsets(ncells_, 0),
        localOffsets(ncells_, 0);
    for (size_t p = 0; p < local; ++p) {
      size_t index[3];
      for (size_t d = 0; d < 3; ++d) {
        double x = std::fmod(buffer.coordinates[3 * (begin + p) + d], dim_);
        if (x < 0.0) x += dim_;
        if (!(x < dim_)) x = 0.0;
        buffer.coordinates[3 * (begin + p) + d] = x;
        index[d] = std::min(static_cast<size_t>(x * iwidth_), cdim_ - 1);
      }
      cellIds[p] = index[2] + cdim_ * (index[1] + cdim_ * index[0]);
      ++localCounts[cellIds[p]];
    }
#ifdef USE_MPI
    if (nranks_ > 1) {
      MPI_Allreduce(localCounts.data(), globalCounts.data(), static_cast<int>(ncells_),
                    MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
      MPI_Exscan(localCounts.data(), preceding.data(), static_cast<int>(ncells_),
                 MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
      if (rank_ == 0) std::fill(preceding.begin(), preceding.end(), 0);
    } else
#endif
    {
      globalCounts = localCounts;
    }
    long long offset = 0;
    size_t localOffset = 0;
    for (size_t c = 0; c < ncells_; ++c) {
      cellOffsets[c] = offset;
      rankOffsets[c] = offset + preceding[c];
      localOffsets[c] = localOffset;
      offset += globalCounts[c];
      localOffset += static_cast<size_t>(localCounts[c]);
    }
    std::vector<size_t> next(localOffsets.begin(), localOffsets.end());
    for (size_t p = 0; p < local; ++p) permutation[next[cellIds[p]]++] = p;

    auto reorder = [&](const std::vector<double> &input, size_t width) {
      std::vector<double> output(local * width);
      for (size_t p = 0; p < local; ++p)
        for (size_t d = 0; d < width; ++d)
          output[p * width + d] = input[(begin + permutation[p]) * width + d];
      return output;
    };
    const std::string group = "/PartType" + std::to_string(type);
    auto coordinates = reorder(buffer.coordinates, 3);
    auto velocities = reorder(buffer.velocities, 3);
    auto masses = reorder(buffer.masses, 1);
    std::vector<unsigned long long> ids(local);
    for (size_t p = 0; p < local; ++p) ids[p] = buffer.ids[begin + permutation[p]];

    auto writeFields = [&](hid_t file, hid_t dxpl) {
      write_cell_ordered_dataset(file, group + "/Coordinates", coordinates, localCounts,
                                 rankOffsets, localOffsets, 3, dxpl);
      write_cell_ordered_dataset(file, group + "/Velocities", velocities, localCounts,
                                 rankOffsets, localOffsets, 3, dxpl);
      write_cell_ordered_dataset(file, group + "/Masses", masses, localCounts,
                                 rankOffsets, localOffsets, 1, dxpl);
      write_cell_ordered_dataset(file, group + "/ParticleIDs", ids, localCounts,
                                 rankOffsets, localOffsets, 1, dxpl);
      if (type == GAS_PARTTYPE && doBaryons) {
        auto smoothing = reorder(buffer.smoothing_lengths, 1);
        auto energy = reorder(buffer.internal_energies, 1);
        write_cell_ordered_dataset(file, group + "/SmoothingLength", smoothing, localCounts,
                                   rankOffsets, localOffsets, 1, dxpl);
        write_cell_ordered_dataset(file, group + "/InternalEnergy", energy, localCounts,
                                   rankOffsets, localOffsets, 1, dxpl);
      }
    };

#ifdef USE_PARALLEL_HDF5
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t file = H5Fopen(fname_.c_str(), H5F_ACC_RDWR, fapl);
    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);
    writeFields(file, dxpl);
    H5Pclose(dxpl);
    H5Fclose(file);
    H5Pclose(fapl);
#else
    for (int writer = 0; writer < nranks_; ++writer) {
#ifdef USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      if (writer != rank_) continue;
      hid_t file = H5Fopen(fname_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
      writeFields(file, H5P_DEFAULT);
      H5Fclose(file);
    }
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

    if (rank_ == 0) {
      const std::string name = "PartType" + std::to_string(type);
      HDFWriteDataset(fname_, "/Cells/Counts/" + name, globalCounts);
      HDFWriteDataset(fname_, "/Cells/OffsetsInFile/" + name, cellOffsets);
      HDFWriteDataset(fname_, "/Cells/Files/" + name, std::vector<int>(ncells_, 0));
      write_unit_attributes(fname_, group + "/Coordinates", 0, 1, 0, 1);
      write_unit_attributes(fname_, group + "/Velocities", 0, 1, -1, 0);
      write_unit_attributes(fname_, group + "/Masses", 1, 0, 0, 0);
      write_unit_attributes(fname_, group + "/ParticleIDs", 0, 0, 0, 0);
      if (type == GAS_PARTTYPE && doBaryons) {
        write_unit_attributes(fname_, group + "/SmoothingLength", 0, 1, 0, 1);
        write_unit_attributes(fname_, group + "/InternalEnergy", 0, 2, -2, 0);
      }
    }
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  /**
   * @brief Have rank zero create the shared file, particle datasets, and
   * /Cells.
   *
   * All particle counts are known from the completed streaming callbacks.
   * Other ranks wait at the finalize barrier before opening the file to write.
   */
  void prepare_output_file() {
    if (rank_ != 0) return;
    HDFCreateFile(fname_);
    for (unsigned type = 0; type < NTYPES; ++type) {
      const particle_buffer &buffer = buffers_[type];
      if (buffer.coordinates.empty()) continue;
      const std::string group = "PartType" + std::to_string(type);
      HDFCreateGroup(fname_, group);
      const size_t count = buffer.coordinates.size() / 3;
      create_particle_dataset<double>("/" + group + "/Coordinates", count, 3);
      create_particle_dataset<double>("/" + group + "/Velocities", count, 3);
      create_particle_dataset<double>("/" + group + "/Masses", count);
      if (useLongIDs)
        create_particle_dataset<unsigned long long>("/" + group + "/ParticleIDs", count);
      else
        create_particle_dataset<unsigned int>("/" + group + "/ParticleIDs", count);
      if (type == GAS_PARTTYPE && doBaryons) {
        create_particle_dataset<double>("/" + group + "/SmoothingLength", count);
        create_particle_dataset<double>("/" + group + "/InternalEnergy", count);
      }
    }
    create_cell_tree(fname_);
  }

  /**
   * @brief Buffer a scalar component from MUSIC's streaming callback.
   *
   * IDs are filled after all particle counts are known; other scalar fields
   * are generated by the mass and gas-property callbacks.
   */
  template <typename T>
  void writeHDF5_a(const std::string &fieldName, int partTypeNum,
                   const std::vector<T> &data) {
    particle_buffer &buffer = buffers_[partTypeNum];
    if (fieldName == "ParticleIDs") {
      buffer.ids.assign(data.begin(), data.end());
    } else if (fieldName == "Masses") {
      buffer.masses.assign(data.begin(), data.end());
    } else if (fieldName == "SmoothingLength") {
      buffer.smoothing_lengths.assign(data.begin(), data.end());
    } else if (fieldName == "InternalEnergy") {
      buffer.internal_energies.assign(data.begin(), data.end());
    } else {
      throw std::runtime_error("Unknown SWIFT scalar field " + fieldName);
    }
  }

  /**
   * @brief Buffer one xyz component, or return a buffered component for mixing.
   *
   * Coarse DM's baryon-relative velocity and displacement corrections read
   * the coarse array before replacing it; readFlag preserves that behavior.
   */
  template <typename T>
  void writeHDF5_b(const std::string &fieldName, int coord, int partTypeNum,
                   std::vector<T> &data, bool readFlag = false) {
    particle_buffer &buffer = buffers_[partTypeNum];
    std::vector<double> &field = fieldName == "Coordinates" ? buffer.coordinates
                                                              : buffer.velocities;
    if (field.empty()) field.resize(3 * data.size());
    if (field.size() != 3 * data.size()) {
      throw std::runtime_error("Inconsistent SWIFT vector field size " + fieldName);
    }
    for (size_t i = 0; i < data.size(); ++i) {
      if (readFlag) data[i] = static_cast<T>(field[3 * i + coord]);
      else field[3 * i + coord] = static_cast<double>(data[i]);
    }
  }

  // called from finalize()
  void generateAndWriteIDs(void) {
    long long offset = 1; // don't use ID==0
    nPartTotAllTypes = 0;

    for (size_t i = 0; i < nPartTotal.size(); i++) {
      if (!nPartTotal[i])
        continue;

      nPartTotAllTypes += nPartTotal[i];

      if (!useLongIDs) {
        std::vector<unsigned int> ids = std::vector<unsigned int>(nPartTotal[i]);
        for (int j = 0; j < nPartTotal[i]; j++)
          ids[j] = offset + j;

        writeHDF5_a("ParticleIDs", i, ids);
      } else {
        std::vector<long long> ids = std::vector<long long>(nPartTotal[i]);
        for (long long j = 0; j < nPartTotal[i]; j++)
          ids[j] = offset + j;

        writeHDF5_a("ParticleIDs", i, ids);
      }

      // make IDs of all particle types sequential (unique) = unnecessary, but consistent with gadget output format
      offset += nPartTotal[i];
    }
  }

  void countLeafCells(const grid_hierarchy &gh) {
    npfine = 0;
    npart = 0;
    npcoarse = 0;

    npfine = gh.count_leaf_cells(gh.levelmax(), gh.levelmax());
    npart = gh.count_leaf_cells(gh.levelmin(), gh.levelmax());

    if (levelmax_ != levelmin_) // multimass
      npcoarse = gh.count_leaf_cells(gh.levelmin(), gh.levelmax() - 1);

    // exact, since the finest level always consists of whole 2x2x2 blocks
    npfine_dm = llround(npfine * pmask.dm_fraction());
    npfine_gas = llround(npfine * pmask.gas_fraction());
  }

  //! Mean gas interparticle spacing (two fine cells for masked loads).
  double gas_smoothing_length() const {
    return boxSize * posFac / std::pow(2, levelmax_) / std::cbrt(pmask.gas_fraction());
  }

  template <typename T> void __write_dm_mass(const grid_hierarchy &gh) {
    countLeafCells(gh);

    // Fine particles: the species' mass in a fine cell, spread over the cells
    // it occupies for masked loads.
    int ilevel = gh.levelmax();
    const double omega_dm = doBaryons ? omega0 - omega_b : omega0;
    const T fine_mass = omega_dm * rhoCrit * pow(boxSize * posFac, 3.0) / pow(2, 3 * levelmax_) /
                        pmask.dm_fraction();

    std::vector<T> data(npfine_dm);
    size_t count = 0;

    for (unsigned i = 0; i < gh.get_grid(ilevel)->size(0); ++i)
      for (unsigned j = 0; j < gh.get_grid(ilevel)->size(1); ++j)
        for (unsigned k = 0; k < gh.get_grid(ilevel)->size(2); ++k)
          if (gh.is_in_mask(ilevel, i, j, k) && !gh.is_refined(ilevel, i, j, k) && pmask.is_dm(i, j, k)) {
            data[count++] = fine_mass;
          }

    writeHDF5_a("Masses", HIGHRES_DM_PARTTYPE, data); // write fine DM

    if (count != npfine_dm)
      throw std::runtime_error("Internal consistency error while writing fine DM masses");

    // fill levelcount for header
    levelcounts = std::vector<size_t>(levelmax_ - levelmin_ + 1);
    for (int ilevel = gh.levelmax(); ilevel >= (int)gh.levelmin(); --ilevel)
      levelcounts[gh.levelmax() - ilevel] = gh.count_leaf_cells(ilevel, ilevel);

    if (levelmax_ != levelmin_) {
      size_t count = 0;
      std::vector<T> data(npcoarse);

      for (int ilevel = gh.levelmax() - 1; ilevel >= (int)gh.levelmin(); --ilevel) {
        // baryon particles live only on finest grid, these particles here are total matter particles
        T pmass = omega0 * rhoCrit * pow(boxSize * posFac, 3.0) / pow(2, 3 * ilevel);

        for (unsigned i = 0; i < gh.get_grid(ilevel)->size(0); ++i)
          for (unsigned j = 0; j < gh.get_grid(ilevel)->size(1); ++j)
            for (unsigned k = 0; k < gh.get_grid(ilevel)->size(2); ++k)
              if (gh.is_in_mask(ilevel, i, j, k) && !gh.is_refined(ilevel, i, j, k)) {
                data[count++] = pmass;
              }
      }

      if (count != npcoarse)
        throw std::runtime_error("Internal consistency error while writing masses");

      writeHDF5_a("Masses", coarsePartType, data); // write DM
    }
  }

  template <typename T> void __write_dm_position(int coord, const grid_hierarchy &gh) {
    countLeafCells(gh);

    // update header
    hsize_t offset_fine = 0, offset_coarse = 0;

    for (unsigned i = 0; i < numFiles; i++) {
      hsize_t dims_fine = ceil(npfine_dm / numFiles);
      hsize_t dims_coarse = ceil(npcoarse / numFiles);

      if (i == numFiles - 1) {
        dims_fine = npfine_dm - offset_fine;
        dims_coarse = npcoarse - offset_coarse;
      }

      nPart[HIGHRES_DM_PARTTYPE] = dims_fine;
      nPart[coarsePartType] = dims_coarse;

      offset_fine += dims_fine;
      offset_coarse += dims_coarse;
    }

    nPartTotal[HIGHRES_DM_PARTTYPE] = npfine_dm;
    nPartTotal[coarsePartType] = npcoarse;

    // FINE: collect displacements and convert to absolute coordinates with correct units
    int ilevel = gh.levelmax();

    std::vector<T> data(npfine_dm);
    size_t count = 0;

    for (unsigned i = 0; i < gh.get_grid(ilevel)->size(0); ++i)
      for (unsigned j = 0; j < gh.get_grid(ilevel)->size(1); ++j)
        for (unsigned k = 0; k < gh.get_grid(ilevel)->size(2); ++k)
          if (gh.is_in_mask(ilevel, i, j, k) && !gh.is_refined(ilevel, i, j, k) && pmask.is_dm(i, j, k)) {
            double xx[3];
            gh.cell_pos(ilevel, i, j, k, xx);

            xx[coord] = (xx[coord] + (*gh.get_grid(ilevel))(i, j, k)) * boxSize;
            xx[coord] = fmod(xx[coord] + boxSize, boxSize); //added fmod here (see below), otherwise only takes latter value

            data[count++] = (T)(xx[coord] * posFac);
          }

    writeHDF5_b("Coordinates", coord, HIGHRES_DM_PARTTYPE, data); // write fine DM

    if (count != npfine_dm)
      throw std::runtime_error("Internal consistency error while writing fine DM pos");

    // COARSE: collect displacements and convert to absolute coordinates with correct units
    if (levelmax_ != levelmin_) // multimass
    {
      data = std::vector<T>(npcoarse, 0.0);
      count = 0;

      for (int ilevel = gh.levelmax() - 1; ilevel >= (int)gh.levelmin(); --ilevel)
        for (unsigned i = 0; i < gh.get_grid(ilevel)->size(0); ++i)
          for (unsigned j = 0; j < gh.get_grid(ilevel)->size(1); ++j)
            for (unsigned k = 0; k < gh.get_grid(ilevel)->size(2); ++k)
              if (gh.is_in_mask(ilevel, i, j, k) && !gh.is_refined(ilevel, i, j, k)) {
                double xx[3];
                gh.cell_pos(ilevel, i, j, k, xx);

                xx[coord] = (xx[coord] + (*gh.get_grid(ilevel))(i, j, k)) * boxSize;

                if (!doBaryons) // if so, we will handle the mod in write_gas_position
                  xx[coord] = fmod(xx[coord] + boxSize, boxSize) * posFac;

                data[count++] = (T)xx[coord];
              }

      if (count != npcoarse)
        throw std::runtime_error("Internal consistency error while writing coarse DM pos");

      writeHDF5_b("Coordinates", coord, coarsePartType, data); // write coarse DM
    }
  }

  template <typename T> void __write_dm_velocity(int coord, const grid_hierarchy &gh) {
    countLeafCells(gh);

    // FINE: collect velocities and convert to correct units
    int ilevel = gh.levelmax();

    std::vector<T> data(npfine_dm);
    size_t count = 0;

    for (unsigned i = 0; i < gh.get_grid(ilevel)->size(0); ++i)
      for (unsigned j = 0; j < gh.get_grid(ilevel)->size(1); ++j)
        for (unsigned k = 0; k < gh.get_grid(ilevel)->size(2); ++k)
          if (gh.is_in_mask(ilevel, i, j, k) && !gh.is_refined(ilevel, i, j, k) && pmask.is_dm(i, j, k)) {
            data[count++] = (T)(*gh.get_grid(ilevel))(i, j, k) * velFac;
          }

    writeHDF5_b("Velocities", coord, HIGHRES_DM_PARTTYPE, data); // write fine DM

    if (count != npfine_dm)
      throw std::runtime_error("Internal consistency error while writing fine DM pos");

    // COARSE: collect velocities and convert to correct units
    if (levelmax_ != levelmin_) // multimass
    {
      data = std::vector<T>(npcoarse, 0.0);
      count = 0;

      for (int ilevel = gh.levelmax() - 1; ilevel >= (int)gh.levelmin(); --ilevel)
        for (unsigned i = 0; i < gh.get_grid(ilevel)->size(0); ++i)
          for (unsigned j = 0; j < gh.get_grid(ilevel)->size(1); ++j)
            for (unsigned k = 0; k < gh.get_grid(ilevel)->size(2); ++k)
              if (gh.is_in_mask(ilevel, i, j, k) && !gh.is_refined(ilevel, i, j, k)) {
                data[count++] = (T)(*gh.get_grid(ilevel))(i, j, k) * velFac;
              }

      if (count != npcoarse)
        throw std::runtime_error("Internal consistency error while writing coarse DM pos");

      writeHDF5_b("Velocities", coord, coarsePartType, data); // write coarse DM
    }
  }

  template <typename T> void __write_gas_velocity(int coord, const grid_hierarchy &gh) {
    countLeafCells(gh);

    std::vector<T> gas_data(npfine_gas + npcoarse); // read/write gas at all levels from the gh
    size_t count = 0;

    for (int ilevel = levelmax_; ilevel >= (int)levelmin_; --ilevel)
      for (unsigned i = 0; i < gh.get_grid(ilevel)->size(0); ++i)
        for (unsigned j = 0; j < gh.get_grid(ilevel)->size(1); ++j)
          for (unsigned k = 0; k < gh.get_grid(ilevel)->size(2); ++k)
            if (gh.is_in_mask(ilevel, i, j, k) && !gh.is_refined(ilevel, i, j, k) &&
                (ilevel != (int)levelmax_ || pmask.is_gas(i, j, k))) {
              gas_data[count++] = (T)(*gh.get_grid(ilevel))(i, j, k) * velFac;
            }

    if (count != npfine_gas + npcoarse)
      throw std::runtime_error("Internal consistency error while writing GAS pos");

    // calculate modified DM velocities if: multimass and baryons present
    if (doBaryons && npcoarse) {
      double facb = omega_b / omega0;
      double facc = (omega0 - omega_b) / omega0;

      std::vector<T> dm_data(npcoarse);

      writeHDF5_b("Velocities", coord, coarsePartType, dm_data, true); // read coarse DM vels

      // overwrite
      for (size_t i = 0; i < npcoarse; i++)
        dm_data[i] = facc * dm_data[i] + facb * gas_data[npfine_gas + i];

      writeHDF5_b("Velocities", coord, coarsePartType, dm_data); // overwrite coarse DM vels
    }                                                            // dm_data deallocated

    // restrict gas_data to fine only and request write
    std::vector<T> data(gas_data.begin() + 0, gas_data.begin() + npfine_gas);

    std::vector<T>().swap(gas_data); // deallocate

    writeHDF5_b("Velocities", coord, GAS_PARTTYPE, data); // write highres gas
  }

  template <typename T> void __write_gas_position(int coord, const grid_hierarchy &gh) {
    countLeafCells(gh);

    // update header (will actually write only gas at levelmax)
    hsize_t offset = 0;

    for (unsigned i = 0; i < numFiles; i++) {
      hsize_t dims = ceil(npfine_gas / numFiles);
      if (i == numFiles - 1)
        dims = npfine_gas - offset;

      nPart[GAS_PARTTYPE] = dims;
      offset += dims;
    }

    nPartTotal[GAS_PARTTYPE] = npfine_gas;

    std::vector<double> gas_data(npfine_gas + npcoarse); // read/write gas at all levels from the gh
    size_t count = 0;

    double h = 1.0 / (1ul << gh.levelmax());

    for (int ilevel = gh.levelmax(); ilevel >= (int)gh.levelmin(); --ilevel)
      for (unsigned i = 0; i < gh.get_grid(ilevel)->size(0); ++i)
        for (unsigned j = 0; j < gh.get_grid(ilevel)->size(1); ++j)
          for (unsigned k = 0; k < gh.get_grid(ilevel)->size(2); ++k)
            if (gh.is_in_mask(ilevel, i, j, k) && !gh.is_refined(ilevel, i, j, k) &&
                (ilevel != (int)gh.levelmax() || pmask.is_gas(i, j, k))) {
              double xx[3];
              gh.cell_pos(ilevel, i, j, k, xx);

              // shift particle positions (this has to be done as the same shift
              // is used when computing the convolution kernel for SPH baryons,
              // on every level); masked loads use the unshifted lattice
              if (!pmask.enabled)
                xx[coord] += 0.5 * h;

              xx[coord] = (xx[coord] + (*gh.get_grid(ilevel))(i, j, k)) * boxSize;

              gas_data[count++] = xx[coord];
            }

    if (count != npfine_gas + npcoarse)
      throw std::runtime_error("Internal consistency error while writing coarse DM pos");

    // calculate modified DM coordinates if: multimass and baryons present
    if (doBaryons && npcoarse) {
      double facb = omega_b / omega0;
      double facc = (omega0 - omega_b) / omega0;

      std::vector<T> dm_data(npcoarse);

      writeHDF5_b("Coordinates", coord, coarsePartType, dm_data, true); // read coarse DM vels

      // overwrite
      for (size_t i = 0; i < npcoarse; i++) {
        dm_data[i] = facc * dm_data[i] + facb * gas_data[npfine_gas + i];
        dm_data[i] = fmod(dm_data[i] + boxSize, boxSize) * posFac;
      }

      writeHDF5_b("Coordinates", coord, coarsePartType, dm_data); // overwrite coarse DM vels
    }

    // restrict gas_data to fine only and request write
    // std::vector<float> data( gas_data.begin() + 0, gas_data.begin() + npfine );

    std::vector<T> data(npfine_gas);

    for (size_t i = 0; i < npfine_gas; i++)
      data[i] = (T)(fmod(gas_data[i] + boxSize, boxSize) * posFac);

    std::vector<double>().swap(gas_data); // deallocate

    writeHDF5_b("Coordinates", coord, GAS_PARTTYPE, data); // write highres gas
  }

  template <typename T> void __write_gas_properties(const grid_hierarchy &gh) {
    countLeafCells(gh);

    std::vector<T> masses(npfine_gas);
    std::vector<T> smoothing_lengths(npfine_gas);
    std::vector<T> internal_energies(npfine_gas);

    T gas_mass = omega_b * rhoCrit * pow(boxSize * posFac, 3.0) / pow(2, 3 * levelmax_) /
                 pmask.gas_fraction();
    T smoothing_length = gas_smoothing_length();

    // calculate internal energy for gas
    double npol = (fabs(1.0 - gamma) > 1e-7) ? 1.0 / (gamma - 1.) : 1.0;
    double astart = 1.0 / (1.0 + redshift);
    double h2 = hubbleParam * hubbleParam;
    double adec = 1.0 / (160.0 * pow(omega_b * h2 / 0.022, 2.0 / 5.0));

    Tini = astart < adec ? Tcmb0 / astart : Tcmb0 / astart / astart * adec;

    const double mu = (Tini > 1.e4) ? 4.0 / (8. - 5. * YHe) : 4.0 / (1. + 3. * (1. - YHe));
    // Use the physical constants assumed internally by SWIFT.
    T internal_energy = 1.380649e-16 / 1.67262192369e-24 * Tini * npol / mu /
                        UnitVelocity_in_cm_per_s / UnitVelocity_in_cm_per_s;

    for (size_t i = 0; i < npfine_gas; i++){
      masses[i] = gas_mass;
      smoothing_lengths[i] = smoothing_length;
      internal_energies[i] = internal_energy;
    }

    writeHDF5_a("Masses", GAS_PARTTYPE, masses); // write gas masses and other data
    writeHDF5_a("SmoothingLength", GAS_PARTTYPE, smoothing_lengths);
    writeHDF5_a("InternalEnergy", GAS_PARTTYPE, internal_energies);
  }

public:
  /**
   * @brief Initialise SWIFT units, cell geometry, and output configuration.
   *
   * `output/top_level_cells` defaults to 32 and must lie in [1, 1290].
   * The upper bound keeps cdim^3 inside the signed 32-bit `nr_cells`
   * attribute. SWIFT output is always one file and double precision;
   * `swift_num_files` and `swift_doubleprec` are retained for configuration
   * compatibility but do not change those requirements.
   */
  swift_output_plugin(config_file &cf) : output_plugin(cf), pmask(cf) {
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks_);
#else
    rank_ = 0;
    nranks_ = 1;
#endif

    // ensure that everyone knows we want to do SPH, implies: bsph=1, bbshift=1, decic_baryons=1
    // -> instead of just writing gas densities (which are here ignored), the gas displacements are also written
    cf.insert_value("setup", "do_SPH", "yes");

    // init header and config parameters
    nPartTotal = std::vector<long long>(NTYPES, 0);
    massTable = std::vector<double>(NTYPES, 0.0);

    coarsePartType = cf.get_value_safe<unsigned>("output", "swift_coarsetype", COARSE_DM_DEFAULT_PARTTYPE);
    UnitLength_in_cm = cf.get_value_safe<double>("output", "swift_unitlength", 3.08567758149e24); // 1.0 Mpc
    UnitMass_in_g = cf.get_value_safe<double>("output", "swift_unitmass", 1.98841e43);             // 1.0e10 solar masses
    UnitVelocity_in_cm_per_s = cf.get_value_safe<double>("output", "swift_unitvel", 1e5);    // 1 km/sec

    omega0 = cf.get_value<double>("cosmology", "Omega_m");
    omega_b = cf.get_value<double>("cosmology", "Omega_b");
    omega_L = cf.get_value<double>("cosmology", "Omega_L");
    redshift = cf.get_value<double>("setup", "zstart");
    boxSize = cf.get_value<double>("setup", "boxlength");
    doBaryons = cf.get_value_safe<bool>("setup", "baryons", false);
    useLongIDs = cf.get_value_safe<bool>("output", "swift_longids", true);
    numFiles = cf.get_value_safe<unsigned>("output", "swift_num_files", 1);
    cdim_ = cf.get_value_safe<size_t>("output", "top_level_cells", 32);
    if (numFiles != 1 && rank_ == 0)
      music::wlog.Print("SWIFT : swift_num_files ignored; one shared file is required by cell metadata");
    numFiles = 1;
    if (cdim_ == 0 || cdim_ > 1290) {
      throw std::runtime_error("top_level_cells must be between 1 and 1290");
    }
    ncells_ = cdim_ * cdim_ * cdim_;

    nPart = std::vector<unsigned long long>(NTYPES, 0);

    // factors which multiply positions and velocities
    time = 1.0 / (1.0 + redshift);
    posFac = 3.08567758149e24 / UnitLength_in_cm; // MUSIC uses Mpc internally, i.e. posFac=1e3 for kpc output
    velFac = boxSize * 1e5 / UnitVelocity_in_cm_per_s;

    // critical density
    rhoCrit = 27.7536609198; // in h^2 1e10 M_sol / Mpc^3, using SWIFT constants
    rhoCrit *= pow(UnitLength_in_cm / 3.08567758149e24, 3.0);
    rhoCrit *= (1.98841e43 / UnitMass_in_g);

    // calculate PMGRID suggestion
    pmgrid = pow(2, levelmin_) * 2; // unigrid
    gridboost = 1;

    if (levelmin_ != levelmax_) {
      vec3_t lxref, x0ref, x1ref;
      double pmgrid_new;

      the_region_generator->get_AABB(x0ref, x1ref, levelmax_); // generalized beyond box
      for (int i = 0; i < 3; i++)
        lxref[i] = x1ref[i] - x0ref[i];

      // fraction box length of the zoom region
      lxref[0] = pow((lxref[0] * lxref[1] * lxref[2]), 0.333);

      pmgrid_new = pow(2, levelmax_) * 2; // to cover entire box at highest resolution
      pmgrid_new *= lxref[0];             // only need to cover a fraction

      if ((gridboost = round(pmgrid_new / pmgrid)) > 1)
        gridboost = pow(2, ceil(log(gridboost) / log(2.0))); // round to nearest, higher power of 2
      if (gridboost == 0)
        gridboost = 1;
    }

    // calculate Tini for gas
    hubbleParam = cf.get_value<double>("cosmology", "H0") / 100.0;
    Tcmb0 = cf.get_value_safe<double>("cosmology", "Tcmb0",
                                     cf.get_value_safe<double>("cosmology", "Tcmb", 2.7255));
    gamma = cf.get_value_safe<double>("cosmology", "gamma", 5.0 / 3.0);
    YHe = cf.get_value_safe<double>("cosmology", "YHe", 0.245421); //from monofonIC Planck2018EE+BAO+SN

    double astart = 1.0 / (1.0 + redshift);
    double h2 = hubbleParam * hubbleParam;
    double adec = 1.0 / (160.0 * pow(omega_b * h2 / 0.022, 2.0 / 5.0));

    rhoCrit *= h2;
    posFac /= hubbleParam;

    // SWIFT cell geometry in the same units as stored coordinates.
    dim_ = boxSize * posFac;
    cell_width_ = dim_ / static_cast<double>(cdim_);
    iwidth_ = static_cast<double>(cdim_) / dim_;

    Tini = astart < adec ? Tcmb0 / astart : Tcmb0 / astart / astart * adec;

    // calculate softening suggestion
    softening = (boxSize * posFac) / pow(2, levelmax_) / 40.0;

    if (coarsePartType == GAS_PARTTYPE || coarsePartType == HIGHRES_DM_PARTTYPE)
      throw std::runtime_error("Error: Specified illegal Swift particle type for coarse particles.");
    if (coarsePartType >= NTYPES)
      throw std::runtime_error("swift_coarsetype must be smaller than 7");
    if (coarsePartType == STAR_PARTTYPE)
      music::wlog.Print("WARNING: Specified coarse particle type will collide with stars if USE_SFR enabled.");

    if (rank_ == 0)
      music::ilog.Print("SWIFT : sorting particles into %zu^3 top-level cells "
                        "of width %.4g output length units",
                        cdim_, cell_width_);
#ifdef USE_PARALLEL_HDF5
    if (rank_ == 0) music::ilog.Print("SWIFT : collective parallel HDF5 output on %d rank(s)", nranks_);
#else
    if (rank_ == 0) music::ilog.Print("SWIFT : serialized single-file HDF5 output on %d rank(s)", nranks_);
#endif
  }

  ~swift_output_plugin() {}

  /* ------------------------------------------------------------------------------- */
  void write_dm_mass(const grid_hierarchy &gh) {
    __write_dm_mass<double>(gh);
  }

  void write_dm_position(int coord, const grid_hierarchy &gh) {
    __write_dm_position<double>(coord, gh);
  }

  void write_dm_velocity(int coord, const grid_hierarchy &gh) {
    __write_dm_velocity<double>(coord, gh);
  }

  void write_dm_density(const grid_hierarchy &gh) { /* skip */
  }

  void write_dm_potential(const grid_hierarchy &gh) { /* skip */
  }

  /* ------------------------------------------------------------------------------- */
  void write_gas_velocity(int coord, const grid_hierarchy &gh) {
    __write_gas_velocity<double>(coord, gh);
  }

  void write_gas_position(int coord, const grid_hierarchy &gh) {
    __write_gas_position<double>(coord, gh);
  }

  void write_gas_properties(const grid_hierarchy &gh) {
    __write_gas_properties<double>(gh);
  }

  void write_gas_density(const grid_hierarchy &gh) {
    // if only saving highres gas, then all gas cells have the same initial mass
    // do not write out densities as we write out displacements
    __write_gas_properties<double>(gh);
  }

  void write_gas_potential(const grid_hierarchy &gh) { /* skip */
  }

  void finalize(void) {
    // generate and add contiguous IDs for each particle type we have written
    generateAndWriteIDs();
    prepare_output_file();
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    for (unsigned type = 0; type < NTYPES; ++type)
      write_particle_type(type);
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (rank_ != 0) return;

    std::vector<unsigned int> nPartTotalLW(nPartTotal.size());
    std::vector<unsigned int> nPartTotalHW(nPartTotal.size());
    for (size_t i = 0; i < nPartTotalHW.size(); i++) {
      nPartTotalHW[i] = (unsigned)((size_t)nPartTotal[i] >> 32);
      nPartTotalLW[i] = (unsigned int)((size_t)nPartTotal[i]);
    }

    // output particle counts
    std::cout << " - Swift : wrote " << nPartTotAllTypes << " particles..." << std::endl;
    for (size_t i = 0; i < nPartTotal.size(); i++)
      std::cout << "    type [" << i << "] : " << std::setw(12) << nPartTotal[i] << std::endl;
    std::cout << std::endl;

    // write final header and provenance (rank zero owns file metadata)
    for (unsigned i = 0; i < numFiles; i++) {
      std::string filename = fname_;
      if (numFiles > 1) {
        std::stringstream s;
        s << "." << i << ".hdf5";
        filename.replace(filename.find(".hdf5"), 5, s.str());

        std::cout << "    " << filename;
        for (size_t j = 0; j < nPart.size(); j++)
          std::cout << " " << std::setw(10) << nPart[j];
        std::cout << std::endl;
      }

      // Write UNITS header using the physical constants assumed internally by SWIFT
      HDFCreateGroup(filename, "Units");
      HDFWriteGroupAttribute(filename, "Units", "Unit mass in cgs (U_M)", std::vector<double>(1, UnitMass_in_g));
      HDFWriteGroupAttribute(filename, "Units", "Unit length in cgs (U_L)", std::vector<double>(1, UnitLength_in_cm));
      HDFWriteGroupAttribute(filename, "Units", "Unit time in cgs (U_t)",
                             std::vector<double>(1, UnitLength_in_cm / UnitVelocity_in_cm_per_s));
      HDFWriteGroupAttribute(filename, "Units", "Unit current in cgs (U_I)", std::vector<double>(1, 1.0));
      HDFWriteGroupAttribute(filename, "Units", "Unit temperature in cgs (U_T)", std::vector<double>(1, 1.0));

      HDFCreateGroup(filename, "Code");
      HDFWriteGroupAttribute(filename, "Code", "Code", std::string("SWIFT"));

      HDFCreateGroup(filename, "Header");

      HDFWriteGroupAttribute(filename, "Header", "NumPart_ThisFile", nPart);
      HDFWriteGroupAttribute(filename, "Header", "NumPart_Total", nPartTotalLW);
      HDFWriteGroupAttribute(filename, "Header", "NumPart_Total_HighWord", nPartTotalHW);
      HDFWriteGroupAttribute(filename, "Header", "NumPartTypes", std::vector<int>(1, NTYPES));
      HDFWriteGroupAttribute(filename, "Header", "MassTable", massTable);
      HDFWriteGroupAttribute(filename, "Header", "InitialMassTable", massTable);
      HDFWriteGroupAttribute(filename, "Header", "BoxSize", std::vector<double>(3, boxSize * posFac));
      HDFWriteGroupAttribute(filename, "Header", "NumFilesPerSnapshot", std::vector<int>(1, numFiles));
      HDFWriteGroupAttribute(filename, "Header", "Dimension", std::vector<int>(1, 3));
      HDFWriteGroupAttribute(filename, "Header", "Time", std::vector<double>(1, time));
      HDFWriteGroupAttribute(filename, "Header", "Redshift", std::vector<double>(1, redshift));
      HDFWriteGroupAttribute(filename, "Header", "Scale-factor", std::vector<double>(1, time));
      HDFWriteGroupAttribute(filename, "Header", "Omega0", omega0); //this might be CDM plus baryons
      HDFWriteGroupAttribute(filename, "Header", "OmegaLambda", omega_L);
      HDFWriteGroupAttribute(filename, "Header", "OmegaBaryon", omega_b);
      HDFWriteGroupAttribute(filename, "Header", "HubbleParam", hubbleParam);
      HDFWriteGroupAttribute(filename, "Header", "Flag_Sfr", 0);
      HDFWriteGroupAttribute(filename, "Header", "Flag_Cooling", 0);
      HDFWriteGroupAttribute(filename, "Header", "Flag_StellarAge", 0);
      HDFWriteGroupAttribute(filename, "Header", "Flag_Metals", 0);
      HDFWriteGroupAttribute(filename, "Header", "Flag_Feedback", 0);
      HDFWriteGroupAttribute(filename, "Header", "Flag_DoublePrecision", 1);
      HDFWriteGroupAttribute(filename, "Header", "Music_levelmin", levelmin_);
      HDFWriteGroupAttribute(filename, "Header", "Music_levelmax", levelmax_);
      HDFWriteGroupAttribute(filename, "Header", "Music_levelcounts", levelcounts);
      HDFWriteGroupAttribute(filename, "Header", "longIDs", (int)useLongIDs);
      HDFWriteGroupAttribute(filename, "Header", "suggested_pmgrid", pmgrid);
      HDFWriteGroupAttribute(filename, "Header", "suggested_highressoft", softening);
      HDFWriteGroupAttribute(filename, "Header", "Flag_Entropy_ICs", std::vector<int>(1, 0));
      HDFWriteGroupAttribute(filename, "Header", "Virtual", 0);
      HDFWriteGroupAttribute(filename, "Header", "ThisFile", i);

      HDFCreateGroup(filename, "ICs_parameters");
      HDFWriteGroupAttribute(filename, "ICs_parameters", "Code", std::string("MUSIC2 - monofonIC"));
      HDFWriteGroupAttribute(filename, "ICs_parameters", "Git Revision", std::string(GIT_REV).empty() ? std::string("unknown") : std::string(GIT_REV));
      HDFWriteGroupAttribute(filename, "ICs_parameters", "Git Tag", std::string(GIT_TAG).empty() ? std::string("unknown") : std::string(GIT_TAG));
      HDFWriteGroupAttribute(filename, "ICs_parameters", "Git Branch", std::string(GIT_BRANCH).empty() ? std::string("unknown") : std::string(GIT_BRANCH));
      HDFWriteGroupAttribute(filename, "ICs_parameters", "Precision", std::string(CMAKE_PRECISION_STR));
      HDFWriteGroupAttribute(filename, "ICs_parameters", "Baryons", (int)doBaryons);
      HDFWriteGroupAttribute(filename, "ICs_parameters", "LPT Order",
                             cf_.get_value_safe<bool>("setup", "use_2LPT", false) ? 2 : 1);
      HDFWriteGroupAttribute(filename, "ICs_parameters", "Transfer Function",
                             cf_.get_value_safe<std::string>("cosmology", "transfer", "unknown"));
      HDFWriteGroupAttribute(filename, "ICs_parameters", "Random Generator",
                             cf_.get_value_safe<std::string>("random", "generator", "MUSIC"));
      HDFWriteGroupAttribute(filename, "ICs_parameters", "Grid Resolution",
                             static_cast<int>(1u << levelmax_));

      write_provenance_metadata(filename);

      HDFCreateGroup(filename, "Cosmology");
      const double unitTime = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
      const double h0Internal = 100.0 * hubbleParam * 1e5 /
                                3.08567758149e24 * unitTime;
      HDFWriteGroupAttribute(filename, "Cosmology", "Cosmological run", std::vector<int>(1, 1));
      HDFWriteGroupAttribute(filename, "Cosmology", "H0 [internal units]", std::vector<double>(1, h0Internal));
      HDFWriteGroupAttribute(filename, "Cosmology", "Omega_b", std::vector<double>(1, omega_b));
      HDFWriteGroupAttribute(filename, "Cosmology", "Omega_cdm", std::vector<double>(1, omega0 - omega_b));
      HDFWriteGroupAttribute(filename, "Cosmology", "Omega_lambda", std::vector<double>(1, omega_L));
      HDFWriteGroupAttribute(filename, "Cosmology", "Omega_m", std::vector<double>(1, omega0));
      const cosmology::parameters &cosmology = the_cosmo_calc->cosmo_param_;
      HDFWriteGroupAttribute(filename, "Cosmology", "Omega_r",
                             std::vector<double>(1, cosmology["Omega_r"]));
      HDFWriteGroupAttribute(filename, "Cosmology", "w_0",
                             std::vector<double>(1, cosmology["w_0"]));
      HDFWriteGroupAttribute(filename, "Cosmology", "w_a",
                             std::vector<double>(1, cosmology["w_a"]));
      HDFWriteGroupAttribute(filename, "Cosmology", "T_CMB_0 [K]",
                             std::vector<double>(1, cosmology["Tcmb"]));
    }

    // give config/parameter file hints
    if (useLongIDs)
      std::cout << " - Swift: Wrote 64bit IDs, enable LONGIDS." << std::endl;
    std::cout << " - Swift: Double precision ICs, set INPUT_IN_DOUBLEPRECISION." << std::endl;
    if (doBaryons)
      std::cout << " - Swift: Wrote high-res gas (only), set REFINEMENT_HIGH_RES_GAS and GENERATE_GAS_IN_ICS with "
                << "SPLIT_PARTICLE_TYPE=" << pow(2, coarsePartType) << "." << std::endl;
    if (levelmax_ != levelmin_)
      std::cout << " - Swift: Have zoom type ICs, set PLACEHIGHRESREGION=" << pow(2, HIGHRES_DM_PARTTYPE)
                << " (suggest PMGRID=" << pmgrid << " with GRIDBOOST=" << gridboost << ")." << std::endl;
    else
      std::cout << " - Swift: Have unigrid type ICs (suggest PMGRID=" << pmgrid << ")." << std::endl;
    if (levelmax_ > levelmin_ + 1)
      std::cout << " - Swift: More than one coarse DM mass using same type, set INDIVIDUAL_GRAVITY_SOFTENING="
                << pow(2, coarsePartType) << " (+" << pow(2, STAR_PARTTYPE) << " if including stars)." << std::endl;
    if (doBaryons)
      std::cout << " - Swift: Set initial gas temperature to " << std::fixed << std::setprecision(3) << Tini << " K."
                << std::endl;
    std::cout << " - Swift: Suggest grav softening = " << std::setprecision(3) << softening << " for high res DM."
              << std::endl;
  }
};

namespace {
output_plugin_creator_concrete<swift_output_plugin> creator("swift");
}

#endif // HAVE_HDF5
