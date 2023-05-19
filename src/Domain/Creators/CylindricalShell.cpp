// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Creators/CylindricalShell.hpp"

#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#include "Domain/BoundaryConditions/Periodic.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/DiscreteRotation.hpp"
#include "Domain/CoordinateMaps/Interval.hpp"
#include "Domain/CoordinateMaps/ProductMaps.hpp"
#include "Domain/CoordinateMaps/ProductMaps.tpp"
#include "Domain/CoordinateMaps/UniformCylindricalEndcap.hpp"
#include "Domain/CoordinateMaps/UniformCylindricalFlatEndcap.hpp"
#include "Domain/CoordinateMaps/UniformCylindricalSide.hpp"
#include "Domain/CoordinateMaps/Wedge.hpp"
#include "Domain/Creators/BinaryCompactObjectHelpers.hpp"
#include "Domain/Creators/ExpandOverBlocks.hpp"
#include "Domain/DomainHelpers.hpp"
#include "Domain/FunctionsOfTime/FixedSpeedCubic.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/QuaternionFunctionOfTime.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Domain/Structure/ExcisionSphere.hpp"
#include "Domain/Structure/OrientationMap.hpp"
#include "NumericalAlgorithms/RootFinding/QuadraticEquation.hpp"

#include <iostream>

namespace {
std::array<double, 3> rotate_to_z_axis(const std::array<double, 3> input) {
  return discrete_rotation(
      OrientationMap<3>{std::array<Direction<3>, 3>{Direction<3>::lower_zeta(),
                                                    Direction<3>::upper_eta(),
                                                    Direction<3>::upper_xi()}},
      input);
}
std::array<double, 3> rotate_from_z_to_x_axis(
    const std::array<double, 3> input) {
  return discrete_rotation(
      OrientationMap<3>{std::array<Direction<3>, 3>{Direction<3>::upper_zeta(),
                                                    Direction<3>::upper_eta(),
                                                    Direction<3>::lower_xi()}},
      input);
}
std::array<double, 3> flip_about_xy_plane(const std::array<double, 3> input) {
  return std::array<double, 3>{input[0], input[1], -input[2]};
}
}  // namespace

namespace domain::creators {
CylindricalShell::CylindricalShell(
    std::array<double, 3> center_A, std::array<double, 3> center_B,
    double radius_A, double radius_B, bool include_inner_sphere_A,
    bool include_inner_sphere_B, double outer_radius, bool use_equiangular_map,
    const typename InitialRefinement::type& initial_refinement,
    const typename InitialGridPoints::type& initial_grid_points,
    std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
        inner_boundary_condition,
    std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
        outer_boundary_condition,
    const Options::Context& context)
    : center_A_(rotate_to_z_axis(center_A)),
      center_B_(rotate_to_z_axis(center_B)),
      radius_A_(radius_A),
      radius_B_(radius_B),
      include_inner_sphere_A_(include_inner_sphere_A),
      outer_radius_(outer_radius),
      use_equiangular_map_(use_equiangular_map),
      inner_boundary_condition_(std::move(inner_boundary_condition)),
      outer_boundary_condition_(std::move(outer_boundary_condition)) {
  std::cout << "Inside CylindricalShell constructor." << std::endl;
  if (center_A_[2] <= 0.0) {
    PARSE_ERROR(
        context,
        "The x-coordinate of the input CenterA is expected to be positive");
  }
  if (radius_A_ <= 0.0) {
    PARSE_ERROR(context, "RadiusA is expected to be positive");
  }

  if ((outer_boundary_condition_ == nullptr) xor
      (inner_boundary_condition_ == nullptr)) {
    PARSE_ERROR(context,
                "Must specify either both inner and outer boundary conditions "
                "or neither.");
  }
  using domain::BoundaryConditions::is_periodic;
  if (is_periodic(inner_boundary_condition_) or
      is_periodic(outer_boundary_condition_)) {
    PARSE_ERROR(
        context,
        "Cannot have periodic boundary conditions with a binary domain");
  }

  // outer_radius_A is the outer radius of the inner sphere A, if it exists.
  // If the inner sphere A does not exist, then outer_radius_A is the same
  // as radius_A_.
  // If the inner sphere does exist, the algorithm for computing
  // outer_radius_A is the same as in SpEC when there is one inner shell.
  outer_radius_A_ = include_inner_sphere_A_
                        ? radius_A_ + 0.5 * (std::abs(center_A_[2]) - radius_A_)
                        : radius_A_;

  number_of_blocks_ = 46;
  if (include_inner_sphere_A) {
    number_of_blocks_ += 14;
  }
  // Create block names and groups
  auto add_filled_cylinder_name = [this](const std::string& prefix,
                                         const std::string& group_name) {
    for (const std::string& where :
         {"Center"s, "East"s, "North"s, "West"s, "South"s}) {
      const std::string name =
          std::string(prefix).append("FilledCylinder").append(where);
      block_names_.push_back(name);
      block_groups_[group_name].insert(name);
    }
  };
  auto add_cylinder_name = [this](const std::string& prefix,
                                  const std::string& group_name) {
    for (const std::string& where : {"East"s, "North"s, "West"s, "South"s}) {
      const std::string name =
          std::string(prefix).append("Cylinder").append(where);
      block_names_.push_back(name);
      block_groups_[group_name].insert(name);
    }
  };

  // EA Filled Cylinder
  // 5 blocks: 9 thru 13
  add_filled_cylinder_name("EA", "InnerA");

  // EA Cylinder
  // 4 blocks: 14 thru 17
  add_cylinder_name("EA", "InnerA");

  // MA Filled Cylinder
  // 5 blocks: 27 thru 31
  add_filled_cylinder_name("MA", "InnerA");

  if (include_inner_sphere_A) {
    // 5 blocks
    add_filled_cylinder_name("InnerSphereEA", "InnerSphereA");
    // 5 blocks
    add_filled_cylinder_name("InnerSphereMA", "InnerSphereA");
    // 4 blocks
    add_cylinder_name("InnerSphereEA", "InnerSphereA");
  }
  // Expand initial refinement over all blocks
  const ExpandOverBlocks<size_t, 3> expand_over_blocks{block_names_,
                                                       block_groups_};
  try {
    initial_refinement_ = std::visit(expand_over_blocks, initial_refinement);
  } catch (const std::exception& error) {
    PARSE_ERROR(context, "Invalid 'InitialRefinement': " << error.what());
  }
  try {
    initial_grid_points_ = std::visit(expand_over_blocks, initial_grid_points);
  } catch (const std::exception& error) {
    PARSE_ERROR(context, "Invalid 'InitialGridPoints': " << error.what());
  }

  // Now we must change the initial refinement and initial grid points
  // for certain blocks, because the [r, theta, perp] directions do
  // not always correspond to [xi, eta, zeta].  The values in
  // initial_refinement_ must correspond to [xi, eta, zeta].
  //
  // In particular, for cylinders: [xi, eta, zeta] = [r, theta, perp]
  // but for filled cylinders: [xi, eta, zeta] = [perp, theta, r].

  auto swap_refinement_and_grid_points_xi_zeta = [this](const size_t block_id) {
    size_t val = gsl::at(initial_refinement_[block_id], 0);
    gsl::at(initial_refinement_[block_id], 0) =
        gsl::at(initial_refinement_[block_id], 2);
    gsl::at(initial_refinement_[block_id], 2) = val;
    val = gsl::at(initial_grid_points_[block_id], 0);
    gsl::at(initial_grid_points_[block_id], 0) =
        gsl::at(initial_grid_points_[block_id], 2);
    gsl::at(initial_grid_points_[block_id], 2) = val;
  };

  // CA Filled Cylinder
  // 5 blocks: 0 thru 4
  for (size_t block = 0; block < 5; ++block) {
    swap_refinement_and_grid_points_xi_zeta(block);
  }

  // Now do the filled cylinders for the inner and outer shells,
  // if they are present.
  size_t current_block = 5;
  if (include_inner_sphere_A) {
    for (size_t block = 0; block < 10; ++block) {
      swap_refinement_and_grid_points_xi_zeta(current_block++);
    }
    current_block += 4;
  }
}

CylindricalShell::CylindricalShell(
    bco::TimeDependentMapOptions time_dependent_options,
    std::array<double, 3> center_A, std::array<double, 3> center_B,
    double radius_A, double radius_B, bool include_inner_sphere_A,
    bool include_inner_sphere_B, double outer_radius, bool use_equiangular_map,
    const typename InitialRefinement::type& initial_refinement,
    const typename InitialGridPoints::type& initial_grid_points,
    std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
        inner_boundary_condition,
    std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
        outer_boundary_condition,
    const Options::Context& context)
    : CylindricalShell(center_A, center_B, radius_A, radius_B,
                       include_inner_sphere_A, include_inner_sphere_B,
                       outer_radius, use_equiangular_map, initial_refinement,
                       initial_grid_points, std::move(inner_boundary_condition),
                       std::move(outer_boundary_condition), context) {
  // The size map, which is applied from the grid to distorted frame, currently
  // needs to start and stop at certain radii around each excision. If the inner
  // spheres aren't included, the outer radii would have to be in the middle of
  // a block. With the inner spheres, the outer radii can be at block
  // boundaries.
  if (not(include_inner_sphere_A)) {
    PARSE_ERROR(context,
                "To use the CylindricalBBH domain with time-dependent maps, "
                "you must include the inner spheres for both objects. "
                "Currently, one or both objects is missing the inner spheres.");
  }

  time_dependent_options_ = std::move(time_dependent_options);

  time_dependent_options_->build_maps(
      std::array{rotate_from_z_to_x_axis(center_A_),
                 rotate_from_z_to_x_axis(center_B_)},
      std::array{std::optional<double>{radius_A_},
                 std::optional<double>{radius_B_}},
      std::array{std::optional<double>{outer_radius_A_},
                 std::optional<double>{outer_radius_B_}},
      outer_radius_);
}

Domain<3> CylindricalShell::create_domain() const {
  std::vector<std::unique_ptr<
      domain::CoordinateMapBase<Frame::BlockLogical, Frame::Inertial, 3>>>
      coordinate_maps{};

  const OrientationMap<3> rotate_to_x_axis{std::array<Direction<3>, 3>{
      Direction<3>::upper_zeta(), Direction<3>::upper_eta(),
      Direction<3>::lower_xi()}};

  const OrientationMap<3> rotate_to_minus_x_axis{std::array<Direction<3>, 3>{
      Direction<3>::lower_zeta(), Direction<3>::upper_eta(),
      Direction<3>::upper_xi()}};

  // The labels EA, EB, EE, etc are from Figure 20 of
  // https://arxiv.org/abs/1206.3015
  //
  // center_EA and radius_EA are the center and outer-radius of the
  // cylindered-sphere EA in Figure 20.
  //
  // center_EB and radius_EB are the center and outer-radius of the
  // cylindered-sphere EB in Figure 20.
  //
  // radius_MB is eq. A16 or A23 in the paper (depending on whether
  // the EE spheres exist), and is the radius of the circle where the EB
  // sphere intersects the cutting plane.
  const std::array<double, 3> center_EA = {0.0, 0.0, center_A_[2]};
  const double radius_MB = std::abs(center_B_[2]);
  const double radius_EA = sqrt(square(center_EA[2]) + square(radius_MB));

  // Construct vector<CoordMap>s that go from logical coordinates to
  // various blocks making up a unit right cylinder.  These blocks are
  // either the central square blocks, or the surrounding wedge
  // blocks. The radii and bounds are what are expected by the
  // UniformCylindricalEndcap maps, (except cylinder_inner_radius, which
  // determines the internal block boundaries inside the cylinder, and
  // which the UniformCylindricalEndcap maps don't care about).
  const double cylinder_inner_radius = 0.5;
  const double cylinder_outer_radius = 1.0;
  const double cylinder_lower_bound_z = -1.0;
  const double cylinder_upper_bound_z = 1.0;
  const auto logical_to_cylinder_center_maps =
      cyl_wedge_coord_map_center_blocks(
          cylinder_inner_radius, cylinder_lower_bound_z, cylinder_upper_bound_z,
          use_equiangular_map_);
  const auto logical_to_cylinder_surrounding_maps =
      cyl_wedge_coord_map_surrounding_blocks(
          cylinder_inner_radius, cylinder_outer_radius, cylinder_lower_bound_z,
          cylinder_upper_bound_z, use_equiangular_map_, 0.0);

  // Lambda that takes a UniformCylindricalEndcap map and a
  // DiscreteRotation map, composes it with the logical-to-cylinder
  // maps, and adds it to the list of coordinate maps. Also adds
  // boundary conditions if requested.
  auto add_endcap_to_list_of_maps =
      [&coordinate_maps, &logical_to_cylinder_center_maps,
       &logical_to_cylinder_surrounding_maps](
          const CoordinateMaps::UniformCylindricalEndcap& endcap_map,
          const CoordinateMaps::DiscreteRotation<3>& rotation_map) {
        auto new_logical_to_cylinder_center_maps =
            domain::make_vector_coordinate_map_base<Frame::BlockLogical,
                                                    Frame::Inertial, 3>(
                logical_to_cylinder_center_maps, endcap_map, rotation_map);
        coordinate_maps.insert(
            coordinate_maps.end(),
            std::make_move_iterator(
                new_logical_to_cylinder_center_maps.begin()),
            std::make_move_iterator(new_logical_to_cylinder_center_maps.end()));
        auto new_logical_to_cylinder_surrounding_maps =
            domain::make_vector_coordinate_map_base<Frame::BlockLogical,
                                                    Frame::Inertial, 3>(
                logical_to_cylinder_surrounding_maps, endcap_map, rotation_map);
        coordinate_maps.insert(
            coordinate_maps.end(),
            std::make_move_iterator(
                new_logical_to_cylinder_surrounding_maps.begin()),
            std::make_move_iterator(
                new_logical_to_cylinder_surrounding_maps.end()));
      };

  // Construct vector<CoordMap>s that go from logical coordinates to
  // various blocks making up a right cylindrical shell of inner radius 1,
  // outer radius 2, and z-extents from -1 to +1.  These blocks are
  // either the central square blocks, or the surrounding wedge
  // blocks. The radii and bounds are what are expected by the
  // UniformCylindricalEndcap maps.
  const double cylindrical_shell_inner_radius = 1.0;
  const double cylindrical_shell_outer_radius = 2.0;
  const double cylindrical_shell_lower_bound_z = -1.0;
  const double cylindrical_shell_upper_bound_z = 1.0;
  const auto logical_to_cylindrical_shell_maps =
      cyl_wedge_coord_map_surrounding_blocks(
          cylindrical_shell_inner_radius, cylindrical_shell_outer_radius,
          cylindrical_shell_lower_bound_z, cylindrical_shell_upper_bound_z,
          use_equiangular_map_, 1.0);

  // Lambda that takes a UniformCylindricalSide map and a DiscreteRotation
  // map, composes it with the logical-to-cylinder maps, and adds it
  // to the list of coordinate maps.  Also adds boundary conditions if
  // requested.
  auto add_side_to_list_of_maps =
      [&coordinate_maps, &logical_to_cylindrical_shell_maps](
          const CoordinateMaps::UniformCylindricalSide& side_map,
          const CoordinateMaps::DiscreteRotation<3>& rotation_map) {
        auto new_logical_to_cylindrical_shell_maps =
            domain::make_vector_coordinate_map_base<Frame::BlockLogical,
                                                    Frame::Inertial, 3>(
                logical_to_cylindrical_shell_maps, side_map, rotation_map);
        coordinate_maps.insert(
            coordinate_maps.end(),
            std::make_move_iterator(
                new_logical_to_cylindrical_shell_maps.begin()),
            std::make_move_iterator(
                new_logical_to_cylindrical_shell_maps.end()));
      };

  // z_cut_CA_lower is the lower z_plane position for the CA endcap,
  // defined by https://arxiv.org/abs/1206.3015 in the bulleted list
  // after Eq. (A.19) EXCEPT that here we use a factor of 1.6 instead of 1.5
  // to put the plane farther from center_A.
  const double z_cut_CA_lower = 1.6 * center_EA[2];
  // z_cut_EA_upper is the upper z_plane position for the EA endcap,
  // which isn't defined in https://arxiv.org/abs/1206.3015 (because the
  // maps are different).  We choose this plane to make the maps
  // less extreme.
  const double z_cut_EA_upper = center_A_[2] + 0.7 * outer_radius_A_;
  // z_cut_EA_lower is the lower z_plane position for the EA endcap,
  // which isn't defined in https://arxiv.org/abs/1206.3015 (because the
  // maps are different).  We choose this plane to make the maps
  // less extreme.
  const double z_cut_EA_lower = center_A_[2] - 0.7 * outer_radius_A_;

  // EA Filled Cylinder
  // 5 blocks: 9 thru 13
  add_endcap_to_list_of_maps(
      CoordinateMaps::UniformCylindricalEndcap(center_A_, center_EA,
                                               outer_radius_A_, radius_EA,
                                               z_cut_EA_upper, z_cut_CA_lower),
      CoordinateMaps::DiscreteRotation<3>(rotate_to_x_axis));

  // EA Cylinder
  // 4 blocks: 14 thru 17
  add_side_to_list_of_maps(
      // For some reason codecov complains about the next line.
      CoordinateMaps::UniformCylindricalSide(  // LCOV_EXCL_LINE
          center_A_, center_EA, outer_radius_A_, radius_EA, z_cut_EA_upper,
          z_cut_EA_lower, z_cut_CA_lower, 0.0),
      CoordinateMaps::DiscreteRotation<3>(rotate_to_x_axis));
  // MA Filled Cylinder
  // 5 blocks: 27 thru 31
  add_endcap_to_list_of_maps(
      CoordinateMaps::UniformCylindricalEndcap(
          flip_about_xy_plane(center_A_), flip_about_xy_plane(center_EA),
          outer_radius_A_, radius_EA, -z_cut_EA_lower, 0.0),
      CoordinateMaps::DiscreteRotation<3>(rotate_to_minus_x_axis));

  if (include_inner_sphere_A_) {
    const double z_cut_upper = center_A_[2] + 0.7 * radius_A_;
    const double z_cut_lower = center_A_[2] - 0.7 * radius_A_;
    // InnerSphereEA Filled Cylinder
    // 5 blocks
    add_endcap_to_list_of_maps(
        // For some reason codecov complains about the next function.
        // LCOV_EXCL_START
        CoordinateMaps::UniformCylindricalEndcap(center_A_, center_A_,
                                                 radius_A_, outer_radius_A_,
                                                 z_cut_upper, z_cut_EA_upper),
        // LCOV_EXCL_START
        CoordinateMaps::DiscreteRotation<3>(rotate_to_x_axis));
    // InnerSphereMA Filled Cylinder
    // 5 blocks
    add_endcap_to_list_of_maps(
        CoordinateMaps::UniformCylindricalEndcap(
            flip_about_xy_plane(center_A_), flip_about_xy_plane(center_A_),
            radius_A_, outer_radius_A_, -z_cut_lower, -z_cut_EA_lower),
        CoordinateMaps::DiscreteRotation<3>(rotate_to_minus_x_axis));
    // InnerSphereEA Cylinder
    // 4 blocks
    add_side_to_list_of_maps(
        // For some reason codecov complains about the next line.
        CoordinateMaps::UniformCylindricalSide(  // LCOV_EXCL_LINE
            center_A_, center_A_, radius_A_, outer_radius_A_, z_cut_upper,
            z_cut_lower, z_cut_EA_upper, z_cut_EA_lower),
        CoordinateMaps::DiscreteRotation<3>(rotate_to_x_axis));
  }
  // Excision spheres
  std::unordered_map<std::string, ExcisionSphere<3>> excision_spheres{};

  std::unordered_map<size_t, Direction<3>> abutting_directions_A;
  size_t first_inner_sphere_block = 0;
  if (include_inner_sphere_A_) {
    for (size_t i = 0; i < 10; ++i) {
      // LCOV_EXCL_START
      abutting_directions_A.emplace(first_inner_sphere_block + i,
                                    Direction<3>::lower_zeta());
      // LCOV_EXCL_STOP
    }
    for (size_t i = 0; i < 4; ++i) {
      // LCOV_EXCL_START
      abutting_directions_A.emplace(first_inner_sphere_block + 10 + i,
                                    Direction<3>::lower_xi());
      // LCOV_EXCL_STOP
    }
    // Block numbers of sphereB might depend on whether there is an inner
    // sphereA layer, so increment here to get that right.
    first_inner_sphere_block += 14;
  } else {
    for (size_t i = 0; i < 5; ++i) {
      abutting_directions_A.emplace(9 + i, Direction<3>::lower_zeta());
      abutting_directions_A.emplace(27 + i, Direction<3>::lower_zeta());
    }
    for (size_t i = 0; i < 4; ++i) {
      abutting_directions_A.emplace(14 + i, Direction<3>::lower_xi());
    }
  }
  excision_spheres.emplace(
      "ExcisionSphereA",
      ExcisionSphere<3>{
          radius_A_,
          tnsr::I<double, 3, Frame::Grid>(rotate_from_z_to_x_axis(center_A_)),
          abutting_directions_A});

  std::unordered_map<size_t, Direction<3>> abutting_directions_B;
  for (size_t i = 0; i < 5; ++i) {
    abutting_directions_B.emplace(18 + i, Direction<3>::lower_zeta());
    abutting_directions_B.emplace(32 + i, Direction<3>::lower_zeta());
    }
    for (size_t i = 0; i < 4; ++i) {
      abutting_directions_B.emplace(23 + i, Direction<3>::lower_xi());
    }
    excision_spheres.emplace(
        "ExcisionSphereB",
        ExcisionSphere<3>{
            radius_B_,
            tnsr::I<double, 3, Frame::Grid>(rotate_from_z_to_x_axis(center_B_)),
            abutting_directions_B});

    Domain<3> domain{std::move(coordinate_maps), std::move(excision_spheres),
                     block_names_, block_groups_};
    return domain;
}

std::vector<DirectionMap<
    3, std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>>>
CylindricalShell::external_boundary_conditions() const {
  if (outer_boundary_condition_ == nullptr) {
    return {};
  }
  std::vector<DirectionMap<
      3, std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>>>
      boundary_conditions{number_of_blocks_};
  /*
    for (size_t i = 0; i < 5; ++i) {
      // if (not include_outer_sphere_) {
      //  CA Filled Cylinder
      boundary_conditions[i][Direction<3>::upper_zeta()] =
          outer_boundary_condition_->get_clone();
      if (not include_inner_sphere_A_) {
        // EA Filled Cylinder
        boundary_conditions[i + 9][Direction<3>::lower_zeta()] =
            inner_boundary_condition_->get_clone();
        // MA Filled Cylinder
        boundary_conditions[i + 27][Direction<3>::lower_zeta()] =
            inner_boundary_condition_->get_clone();
      }
    }
    for (size_t i = 0; i < 4; ++i) {
      //  CA Cylinder
      boundary_conditions[i + 5][Direction<3>::upper_xi()] =
          outer_boundary_condition_->get_clone();
    }

    size_t last_block = 46;
    if (include_inner_sphere_A_) {
      for (size_t i = 0; i < 5; ++i) {
        // InnerSphereEA Filled Cylinder
        boundary_conditions[last_block + i][Direction<3>::lower_zeta()] =
            inner_boundary_condition_->get_clone();
        // InnerSphereMA Filled Cylinder
        boundary_conditions[last_block + i + 5][Direction<3>::lower_zeta()] =
            inner_boundary_condition_->get_clone();
      }
      for (size_t i = 0; i < 4; ++i) {
        // InnerSphereEA Cylinder
        boundary_conditions[last_block + i + 10][Direction<3>::lower_xi()] =
            inner_boundary_condition_->get_clone();
      }
      last_block += 14;
    }
  */
  return boundary_conditions;
}

std::vector<std::array<size_t, 3>> CylindricalShell::initial_extents() const {
  return initial_grid_points_;
}

std::vector<std::array<size_t, 3>> CylindricalShell::initial_refinement_levels()
    const {
  return initial_refinement_;
}

std::unordered_map<std::string,
                   std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
CylindricalShell::functions_of_time(
    const std::unordered_map<std::string, double>& initial_expiration_times)
    const {
  return time_dependent_options_.has_value()
             ? time_dependent_options_->create_functions_of_time(
                   initial_expiration_times)
             : std::unordered_map<
                   std::string,
                   std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>{};
}
}  // namespace domain::creators
