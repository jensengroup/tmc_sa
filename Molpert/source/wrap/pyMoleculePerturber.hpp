#ifndef _PY_MOLECULE_PERTURBER_HPP_
#define _PY_MOLECULE_PERTURBER_HPP_

#include "MoleculePerturber.hpp"
#include "pyMolecularPerturbations.hpp"
#include "pySTL.hpp"
#include <boost/python.hpp>

namespace python = boost::python;


template <class T>
python::object PythonOptional(const std::shared_ptr<T>& ptr) {
  if (!ptr) {
    return python::object(); // Converts to None in Python
  };
  return python::object(*ptr);
};


std::shared_ptr<std::mt19937> PRNGFactory() {
  std::random_device rd;
  return std::shared_ptr<std::mt19937>(new std::mt19937(rd()));
};

void SeedPRNG(std::mt19937& prng, std::uint32_t seed) {
  prng.seed(seed);
};


python::object QueueFront(const MolecularPerturbationQueue& queue) {
  return PythonPerturbation(queue.front().get());
};


python::object ChangeAtomicNumber(
  const MoleculePerturber& perturber,
  const RDKit::ROMol& molecule,
  AtomIdx atom_idx,
  std::mt19937& prng,
  const MolecularConstraints* constraints = nullptr) {
  return PythonOptional(
    perturber.ChangeAtomicNumber(molecule, atom_idx, prng, constraints));
};

python::object ChangeAtomicNumber(
  const MoleculePerturber& perturber,
  const RDKit::ROMol& molecule,
  std::mt19937& prng,
  const MolecularConstraints* constraints = nullptr) {
  return PythonOptional(
    perturber.ChangeAtomicNumber(molecule, prng, constraints));
};

void AtomicNumberChanges(
  const MoleculePerturber& perturber,
  MolecularPerturbationQueue& queue,
  const RDKit::ROMol& molecule,
  AtomIdx atom_idx,
  const MolecularConstraints* constraints = nullptr,
  const python::object& allowed_values = python::object()) {
  // Everything is allowed
  if (allowed_values.is_none()) {
    perturber.AtomicNumberChanges(queue, molecule, atom_idx, constraints);
    return;
  };
  auto values = to_vector<std::uint8_t>(allowed_values);
  perturber.AtomicNumberChanges(queue, molecule, atom_idx, constraints, &values);
};

void AtomicNumberChanges(
  const MoleculePerturber& perturber,
  MolecularPerturbationQueue& queue,
  const RDKit::ROMol& molecule,
  const MolecularConstraints* constraints = nullptr,
  const python::object& allowed_values = python::object()) {
  if (allowed_values.is_none()) {
    perturber.AtomicNumberChanges(queue, molecule, constraints);
    return;
  };
  auto values = to_vector<std::uint8_t>(allowed_values);
  perturber.AtomicNumberChanges(queue, molecule, constraints, &values);
};


python::object ChangeFormalCharge(
  const MoleculePerturber& perturber,
  const RDKit::ROMol& molecule,
  AtomIdx atom_idx,
  std::mt19937& prng,
  const MolecularConstraints* constraints = nullptr) {
  return PythonOptional(
    perturber.ChangeFormalCharge(molecule, atom_idx, prng, constraints));
};

python::object ChangeFormalCharge(
  const MoleculePerturber& perturber,
  const RDKit::ROMol& molecule,
  std::mt19937& prng,
  const MolecularConstraints* constraints = nullptr) {
  return PythonOptional(
    perturber.ChangeFormalCharge(molecule, prng, constraints));
};

void FormalChargeChanges(
  const MoleculePerturber& perturber,
  MolecularPerturbationQueue& queue,
  const RDKit::ROMol& molecule,
  AtomIdx atom_idx,
  const MolecularConstraints* constraints = nullptr,
  const python::object& allowed_values = python::object()) {
  if (allowed_values.is_none()) {
    perturber.FormalChargeChanges(queue, molecule, atom_idx, constraints);
    return;
  };
  auto values = to_vector<std::int8_t>(allowed_values);
  perturber.FormalChargeChanges(queue, molecule, atom_idx, constraints, &values);
};

void FormalChargeChanges(
  const MoleculePerturber& perturber,
  MolecularPerturbationQueue& queue,
  const RDKit::ROMol& molecule,
  const MolecularConstraints* constraints = nullptr,
  const python::object& allowed_values = python::object()) {
  if (allowed_values.is_none()) {
    perturber.FormalChargeChanges(queue, molecule, constraints);
    return;
  };
  auto values = to_vector<std::int8_t>(allowed_values);
  perturber.FormalChargeChanges(queue, molecule, constraints, &values);
};


python::object ChangeExplicitHydrogens(
  const MoleculePerturber& perturber,
  const RDKit::ROMol& molecule,
  AtomIdx atom_idx,
  std::mt19937& prng,
  const MolecularConstraints* constraints = nullptr) {
  return PythonOptional(
    perturber.ChangeExplicitHydrogens(molecule, atom_idx, prng, constraints));
};

python::object ChangeExplicitHydrogens(
  const MoleculePerturber& perturber,
  const RDKit::ROMol& molecule,
  std::mt19937& prng,
  const MolecularConstraints* constraints = nullptr) {
  return PythonOptional(
    perturber.ChangeExplicitHydrogens(molecule, prng, constraints));
};

void ExplicitHydrogenChanges(
  const MoleculePerturber& perturber,
  MolecularPerturbationQueue& queue,
  const RDKit::ROMol& molecule,
  AtomIdx atom_idx,
  const MolecularConstraints* constraints = nullptr,
  const python::object& allowed_values = python::object()) {
  if (allowed_values.is_none()) {
    perturber.ExplicitHydrogenChanges(queue, molecule, atom_idx, constraints);
    return;
  };
  auto values = to_vector<std::uint8_t>(allowed_values);
  perturber.ExplicitHydrogenChanges(queue, molecule, atom_idx, constraints, &values);
};

void ExplicitHydrogenChanges(
  const MoleculePerturber& perturber,
  MolecularPerturbationQueue& queue,
  const RDKit::ROMol& molecule,
  const MolecularConstraints* constraints = nullptr,
  const python::object& allowed_values = python::object()) {
  if (allowed_values.is_none()) {
    perturber.ExplicitHydrogenChanges(queue, molecule, constraints);
    return;
  };
  auto values = to_vector<std::uint8_t>(allowed_values);
  perturber.ExplicitHydrogenChanges(queue, molecule, constraints, &values);
};


python::object ChangeBondType(
  const MoleculePerturber& perturber,
  const RDKit::ROMol& molecule,
  BondIdx bond_idx,
  std::mt19937& prng,
  const MolecularConstraints* constraints = nullptr) {
  return PythonOptional(
    perturber.ChangeBondType(molecule, bond_idx, prng, constraints));
};

python::object ChangeBondType(
  const MoleculePerturber& perturber,
  const RDKit::ROMol& molecule,
  std::mt19937& prng,
  const MolecularConstraints* constraints = nullptr) {
  return PythonOptional(
    perturber.ChangeBondType(molecule, prng, constraints));
};

void BondTypeChanges(
  const MoleculePerturber& perturber,
  MolecularPerturbationQueue& queue,
  const RDKit::ROMol& molecule,
  BondIdx bond_idx,
  const MolecularConstraints* constraints = nullptr,
  const python::object& allowed_values = python::object()) {
  if (allowed_values.is_none()) {
    perturber.BondTypeChanges(queue, molecule, bond_idx, constraints);
    return;
  };
  auto values = to_vector<RDKit::Bond::BondType>(allowed_values);
  perturber.BondTypeChanges(queue, molecule, bond_idx, constraints, &values);
};

void BondTypeChanges(
  const MoleculePerturber& perturber,
  MolecularPerturbationQueue& queue,
  const RDKit::ROMol& molecule,
  const MolecularConstraints* constraints = nullptr,
  const python::object& allowed_values = python::object()) {
  if (allowed_values.is_none()) {
    perturber.BondTypeChanges(queue, molecule, constraints);
    return;
  };
  auto values = to_vector<RDKit::Bond::BondType>(allowed_values);
  perturber.BondTypeChanges(queue, molecule, constraints, &values);
};


python::object InsertAtom(
  const MoleculePerturber& perturber,
  const RDKit::ROMol& molecule,
  const python::object& neighbor_atom_indices,
  std::mt19937& prng,
  const MolecularConstraints* constraints = nullptr) {
  return PythonOptional(
    perturber.InsertAtom(
      molecule, to_vector<AtomIdx>(neighbor_atom_indices), prng, constraints));
};

python::object InsertAtom(
  const MoleculePerturber& perturber,
  const RDKit::ROMol& molecule,
  AtomIdx central_atom_idx,
  std::mt19937& prng,
  const MolecularConstraints* constraints = nullptr) {
  return PythonOptional(
    perturber.InsertAtom(molecule, central_atom_idx, prng, constraints));
};

python::object InsertAtom(
  const MoleculePerturber& perturber,
  const RDKit::ROMol& molecule,
  std::mt19937& prng,
  const MolecularConstraints* constraints = nullptr) {
  return PythonOptional(
    perturber.InsertAtom(molecule, prng, constraints));
};

void AtomInsertions(
  const MoleculePerturber& perturber,
  MolecularPerturbationQueue& queue,
  const RDKit::ROMol& molecule,
  const python::object& neighbor_atom_indices,
  const MolecularConstraints* constraints = nullptr) {
  perturber.AtomInsertions(
    queue, molecule, to_vector<AtomIdx>(neighbor_atom_indices), constraints);
};


python::object DeleteAtom(
  const MoleculePerturber& perturber,
  const RDKit::ROMol& molecule,
  AtomIdx atom_idx,
  std::mt19937& prng,
  const MolecularConstraints* constraints = nullptr) {
  return PythonOptional(
    perturber.DeleteAtom(molecule, atom_idx, prng, constraints));
};

python::object DeleteAtom(
  const MoleculePerturber& perturber,
  const RDKit::ROMol& molecule,
  std::mt19937& prng,
  const MolecularConstraints* constraints = nullptr) {
  return PythonOptional(
    perturber.DeleteAtom(molecule, prng, constraints));
};


python::object InsertBond(
  const MoleculePerturber& perturber,
  const RDKit::ROMol& molecule,
  AtomIdx begin_atom_idx,
  std::mt19937& prng,
  const MolecularConstraints* constraints = nullptr) {
  return PythonOptional(
    perturber.InsertBond(molecule, begin_atom_idx, prng, constraints));
};

python::object InsertBond(
  const MoleculePerturber& perturber,
  const RDKit::ROMol& molecule,
  std::mt19937& prng,
  const MolecularConstraints* constraints = nullptr) {
  return PythonOptional(
    perturber.InsertBond(molecule, prng, constraints));
};


python::object DeleteBond(
  const MoleculePerturber& perturber,
  const RDKit::ROMol& molecule,
  BondIdx bond_idx,
  std::mt19937& prng,
  const MolecularConstraints* constraints = nullptr) {
  return PythonOptional(
    perturber.DeleteBond(molecule, bond_idx, prng, constraints));
};

python::object DeleteBond(
  const MoleculePerturber& perturber,
  const RDKit::ROMol& molecule,
  std::mt19937& prng,
  const MolecularConstraints* constraints = nullptr) {
  return PythonOptional(
    perturber.DeleteBond(molecule, prng, constraints));
};


python::object Perturbation(
  const MoleculePerturber& perturber,
  const RDKit::ROMol& molecule,
  std::mt19937& prng,
  const MolecularConstraints* constraints = nullptr) {
  return PythonPerturbation(perturber(molecule, prng, constraints).get());
};

void Perturbations(
  const MoleculePerturber& perturber,
  MolecularPerturbationQueue& queue,
  const RDKit::ROMol& molecule,
  const MolecularConstraints* constraints = nullptr) {
  perturber(queue, molecule, constraints);
};


template <class T>
python::tuple GetValuesAndWeights(
  const std::vector<T>& values,
  const std::vector<double>& weights) {
  return python::make_tuple(to_list(values), to_list(weights));
};

python::tuple GetAtomicNumbers(
  const MoleculePerturber& perturber) {
  return GetValuesAndWeights(
    perturber.atomic_numbers, perturber.atomic_numbers_weights);
};

python::tuple GetRingAtomicNumbers(
  const MoleculePerturber& perturber) {
  return GetValuesAndWeights(
    perturber.ring_atomic_numbers, perturber.ring_atomic_numbers_weights);
};

python::tuple GetFormalCharges(
  const MoleculePerturber& perturber) {
  return GetValuesAndWeights(
    perturber.formal_charges, perturber.formal_charges_weights);
};

python::tuple GetNExplicitHydrogens(
  const MoleculePerturber& perturber) {
  return GetValuesAndWeights(
    perturber.n_explicit_hydrogens, perturber.n_explicit_hydrogens_weights);
};

python::tuple GetBondTypes(
  const MoleculePerturber& perturber) {
  return GetValuesAndWeights(
    perturber.bond_types, perturber.bond_types_weights);
};

python::tuple GetRingBondTypes(
  const MoleculePerturber& perturber) {
  return GetValuesAndWeights(
    perturber.ring_bond_types, perturber.ring_bond_types_weights);
};

python::tuple GetPerturbationTypes(
  const MoleculePerturber& perturber) {
  python::list types;
  python::list weights;
  for (std::size_t t = 0; t < MolecularPerturbation::n_types; ++t) {
    if (perturber.perturbation_types[t]) {
      types.append(MolecularPerturbation::Type(t));
      weights.append(perturber.perturbation_types_weights[t]);
    };
  };
  return python::make_tuple(types, weights);
};


template <class T>
void SetValuesAndWeights(
  std::vector<T>& values,
  std::vector<double>& weights,
  const python::list& new_values,
  const python::list& new_weights) {
  if (python::len(new_values) != python::len(new_weights)) {
    throw std::length_error("Size mismatch between values and weights");
  };
  values = to_vector<T>(new_values);
  weights = to_vector<double>(new_weights);
};

void SetPerturbationTypes(
  MoleculePerturber& perturber,
  const python::list& perturbation_types,
  const python::list& weights,
  bool reset = true) {
  if (reset) {
    perturber.perturbation_types.reset();
    std::fill(
      perturber.perturbation_types_weights.begin(),
      perturber.perturbation_types_weights.end(), 0.0);
  };
  std::size_t n = python::len(perturbation_types);
  if (python::len(weights) != n) {
    throw std::length_error("Size mismatch between values and weights");
  };
  for (std::size_t i = 0; i < n; ++i) {
    MolecularPerturbation::Type type =
      python::extract<MolecularPerturbation::Type>(perturbation_types[i]);
    double weight = python::extract<double>(weights[i]);
    perturber.perturbation_types.set(type);
    perturber.perturbation_types_weights[type] = weight;
  };
};

void SetAtomicNumbers(
  MoleculePerturber& perturber,
  const python::list& atomic_numbers,
  const python::list& weights) {
  SetValuesAndWeights(
    perturber.atomic_numbers, perturber.atomic_numbers_weights,
    atomic_numbers, weights);
};

void SetRingAtomicNumbers(
  MoleculePerturber& perturber,
  const python::list& atomic_numbers,
  const python::list& weights) {
  SetValuesAndWeights(
    perturber.ring_atomic_numbers, perturber.ring_atomic_numbers_weights,
    atomic_numbers, weights);
};

void SetFormalCharges(
  MoleculePerturber& perturber,
  const python::list& formal_charges,
  const python::list& weights) {
  SetValuesAndWeights(
    perturber.formal_charges, perturber.formal_charges_weights,
    formal_charges, weights);
};

void SetNExplicitHydrogens(
  MoleculePerturber& perturber,
  const python::list& n_explicit_hydrogens,
  const python::list& weights) {
  SetValuesAndWeights(
    perturber.n_explicit_hydrogens, perturber.n_explicit_hydrogens_weights,
    n_explicit_hydrogens, weights);
};

void SetBondTypes(
  MoleculePerturber& perturber,
  const python::list& bond_types,
  const python::list& weights) {
  SetValuesAndWeights(
    perturber.bond_types, perturber.bond_types_weights,
    bond_types, weights);
};

void SetRingBondTypes(
  MoleculePerturber& perturber,
  const python::list& bond_types,
  const python::list& weights) {
  SetValuesAndWeights(
    perturber.ring_bond_types, perturber.ring_bond_types_weights,
    bond_types, weights);
};


void WrapMoleculePerturber() {

  python::class_<std::mt19937, std::shared_ptr<std::mt19937>>(
    "PRNG", python::init<std::uint32_t>((python::arg("seed"))))
    .def("__init__", python::make_constructor(&PRNGFactory))
    .def("seed", SeedPRNG, (python::arg("seed")))
    .def("__call__", &std::mt19937::operator());

  python::class_<MolecularPerturbationQueue>(
    "MolecularPerturbationQueue", python::init<>())
    .def("push", &MolecularPerturbationQueue::push, (
      python::arg("perturbation")))
    .def("pop", &MolecularPerturbationQueue::pop)
    .def("shuffle", &MolecularPerturbationQueue::shuffle)
    .def("clear", &MolecularPerturbationQueue::clear)
    .def("front", QueueFront)
    .def("empty", &MolecularPerturbationQueue::empty)
    .def("__len__", &MolecularPerturbationQueue::size);

  boost::python::pointer_wrapper<const MolecularConstraints*>
    null_constraints (nullptr);

  python::class_<MoleculePerturber>(
    "MoleculePerturber",
    python::init<bool, bool, bool, bool, bool>((
      python::arg("use_chembl_distribution") = true,
      python::arg("use_aromatic_bonds") = false,
      python::arg("acyclic_can_be_aromatic") = false,
      python::arg("cyclicity_based_atomic_numbers") = false,
      python::arg("cyclicity_based_bond_types") = false)))

    .def("SetChEMBLProperties", &MoleculePerturber::SetChEMBLProperties)
    .def("SetAromaticProperties", &MoleculePerturber::SetAromaticProperties, (
      python::arg("acyclic_can_be_aromatic") = false,
      python::arg("use_chembl_distribution") = true))
    .def("SetDecorationSettings", &MoleculePerturber::SetDecorationSettings)
    .def("SetConservativeSettings", &MoleculePerturber::SetConservativeSettings)
    .def("SetBalancedSettings", &MoleculePerturber::SetBalancedSettings)
    .def("SetAdventurousSettings", &MoleculePerturber::SetAdventurousSettings)

    .def<python::object (const MoleculePerturber&, const RDKit::ROMol&, AtomIdx, std::mt19937&, const MolecularConstraints*)>(
      "ChangeAtomicNumber", ChangeAtomicNumber, (
      python::arg("molecule"),
      python::arg("atom_idx"),
      python::arg("prng"),
      python::arg("constraints") = null_constraints))
    .def<python::object (const MoleculePerturber&, const RDKit::ROMol&, std::mt19937&, const MolecularConstraints*)>(
      "ChangeAtomicNumber", ChangeAtomicNumber, (
      python::arg("molecule"),
      python::arg("prng"),
      python::arg("constraints") = null_constraints))
    .def<void (const MoleculePerturber&, MolecularPerturbationQueue&, const RDKit::ROMol&, AtomIdx, const MolecularConstraints*, const python::object&)>(
      "AtomicNumberChanges", &AtomicNumberChanges, (
      python::arg("queue"),
      python::arg("molecule"),
      python::arg("atom_idx"),
      python::arg("constraints") = null_constraints,
      python::arg("allowed_values") = python::object()))
    .def<void (const MoleculePerturber&, MolecularPerturbationQueue&, const RDKit::ROMol&, const MolecularConstraints*, const python::object&)>(
      "AtomicNumberChanges", &AtomicNumberChanges, (
      python::arg("queue"),
      python::arg("molecule"),
      python::arg("constraints") = null_constraints,
      python::arg("allowed_values") = python::object()))

    .def<python::object (const MoleculePerturber&, const RDKit::ROMol&, AtomIdx, std::mt19937&, const MolecularConstraints*)>(
      "ChangeFormalCharge", ChangeFormalCharge, (
      python::arg("molecule"),
      python::arg("atom_idx"),
      python::arg("prng"),
      python::arg("constraints") = null_constraints))
    .def<python::object (const MoleculePerturber&, const RDKit::ROMol&, std::mt19937&, const MolecularConstraints*)>(
      "ChangeFormalCharge", ChangeFormalCharge, (
      python::arg("molecule"),
      python::arg("prng"),
      python::arg("constraints") = null_constraints))
    .def<void (const MoleculePerturber&, MolecularPerturbationQueue&, const RDKit::ROMol&, AtomIdx, const MolecularConstraints*, const python::object&)>(
      "FormalChargeChanges", &FormalChargeChanges, (
      python::arg("queue"),
      python::arg("molecule"),
      python::arg("atom_idx"),
      python::arg("constraints") = null_constraints,
      python::arg("allowed_values") = python::object()))
    .def<void (const MoleculePerturber&, MolecularPerturbationQueue&, const RDKit::ROMol&, const MolecularConstraints*, const python::object&)>(
      "FormalChargeChanges", &FormalChargeChanges, (
      python::arg("queue"),
      python::arg("molecule"),
      python::arg("constraints") = null_constraints,
      python::arg("allowed_values") = python::object()))

    .def<python::object (const MoleculePerturber&, const RDKit::ROMol&, AtomIdx, std::mt19937&, const MolecularConstraints*)>(
      "ChangeExplicitHydrogens", ChangeExplicitHydrogens, (
      python::arg("molecule"),
      python::arg("atom_idx"),
      python::arg("prng"),
      python::arg("constraints") = null_constraints))
    .def<python::object (const MoleculePerturber&, const RDKit::ROMol&, std::mt19937&, const MolecularConstraints*)>(
      "ChangeExplicitHydrogens", ChangeExplicitHydrogens, (
      python::arg("molecule"),
      python::arg("prng"),
      python::arg("constraints") = null_constraints))
    .def<void (const MoleculePerturber&, MolecularPerturbationQueue&, const RDKit::ROMol&, AtomIdx, const MolecularConstraints*, const python::object&)>(
      "ExplicitHydrogenChanges", &ExplicitHydrogenChanges, (
      python::arg("queue"),
      python::arg("molecule"),
      python::arg("atom_idx"),
      python::arg("constraints") = null_constraints,
      python::arg("allowed_values") = python::object()))
    .def<void (const MoleculePerturber&, MolecularPerturbationQueue&, const RDKit::ROMol&, const MolecularConstraints*, const python::object&)>(
      "ExplicitHydrogenChanges", &ExplicitHydrogenChanges, (
      python::arg("queue"),
      python::arg("molecule"),
      python::arg("constraints") = null_constraints,
      python::arg("allowed_values") = python::object()))

    .def<python::object (const MoleculePerturber&, const RDKit::ROMol&, BondIdx, std::mt19937&, const MolecularConstraints*)>(
      "ChangeBondType", ChangeBondType, (
      python::arg("molecule"),
      python::arg("bond_idx"),
      python::arg("prng"),
      python::arg("constraints") = null_constraints))
    .def<python::object (const MoleculePerturber&, const RDKit::ROMol&, std::mt19937&, const MolecularConstraints*)>(
      "ChangeBondType", ChangeBondType, (
      python::arg("molecule"),
      python::arg("prng"),
      python::arg("constraints") = null_constraints))
    .def<void (const MoleculePerturber&, MolecularPerturbationQueue&, const RDKit::ROMol&, BondIdx, const MolecularConstraints*, const python::object&)>(
      "BondTypeChanges", &BondTypeChanges, (
      python::arg("queue"),
      python::arg("molecule"),
      python::arg("bond_idx"),
      python::arg("constraints") = null_constraints,
      python::arg("allowed_values") = python::object()))
    .def<void (const MoleculePerturber&, MolecularPerturbationQueue&, const RDKit::ROMol&, const MolecularConstraints*, const python::object&)>(
      "BondTypeChanges", &BondTypeChanges, (
      python::arg("queue"),
      python::arg("molecule"),
      python::arg("constraints") = null_constraints,
      python::arg("allowed_values") = python::object()))

    .def<python::object (const MoleculePerturber&, const RDKit::ROMol&, const python::object&, std::mt19937&, const MolecularConstraints*)>(
      "InsertAtom", InsertAtom, (
      python::arg("molecule"),
      python::arg("neighbor_atom_indices"),
      python::arg("prng"),
      python::arg("constraints") = null_constraints))
    .def<python::object (const MoleculePerturber&, const RDKit::ROMol&, AtomIdx, std::mt19937&, const MolecularConstraints*)>(
      "InsertAtom", InsertAtom, (
      python::arg("molecule"),
      python::arg("central_atom_idx"),
      python::arg("prng"),
      python::arg("constraints") = null_constraints))
    .def<python::object (const MoleculePerturber&, const RDKit::ROMol&, std::mt19937&, const MolecularConstraints*)>(
      "InsertAtom", InsertAtom, (
      python::arg("molecule"),
      python::arg("prng"),
      python::arg("constraints") = null_constraints))
    .def<void (const MoleculePerturber&, MolecularPerturbationQueue&, const RDKit::ROMol&, const python::object&, const MolecularConstraints*)>(
      "AtomInsertions", AtomInsertions, (
      python::arg("queue"),
      python::arg("molecule"),
      python::arg("neighbor_atom_indices"),
      python::arg("constraints") = null_constraints))
    .def<void (MoleculePerturber::*)(MolecularPerturbationQueue&, const RDKit::ROMol&, AtomIdx, const MolecularConstraints*) const>(
      "AtomInsertions", &MoleculePerturber::AtomInsertions, (
      python::arg("queue"),
      python::arg("molecule"),
      python::arg("central_atom_idx"),
      python::arg("constraints") = null_constraints))
    .def<void (MoleculePerturber::*)(MolecularPerturbationQueue&, const RDKit::ROMol&, const MolecularConstraints*) const>(
      "AtomInsertions", &MoleculePerturber::AtomInsertions, (
      python::arg("queue"),
      python::arg("molecule"),
      python::arg("constraints") = null_constraints))

    .def<python::object (const MoleculePerturber&, const RDKit::ROMol&, AtomIdx, std::mt19937&, const MolecularConstraints*)>(
      "DeleteAtom", DeleteAtom, (
      python::arg("molecule"),
      python::arg("atom_idx"),
      python::arg("prng"),
      python::arg("constraints") = null_constraints))
    .def<python::object (const MoleculePerturber&, const RDKit::ROMol&, std::mt19937&, const MolecularConstraints*)>(
      "DeleteAtom", DeleteAtom, (
      python::arg("molecule"),
      python::arg("prng"),
      python::arg("constraints") = null_constraints))
    .def<void (MoleculePerturber::*)(MolecularPerturbationQueue&, const RDKit::ROMol&, AtomIdx, const MolecularConstraints*) const>(
      "AtomDeletions", &MoleculePerturber::AtomDeletions, (
      python::arg("queue"),
      python::arg("molecule"),
      python::arg("atom_idx"),
      python::arg("constraints") = null_constraints))
    .def<void (MoleculePerturber::*)(MolecularPerturbationQueue&, const RDKit::ROMol&, const MolecularConstraints*) const>(
      "AtomDeletions", &MoleculePerturber::AtomDeletions, (
      python::arg("queue"),
      python::arg("molecule"),
      python::arg("constraints") = null_constraints))

    .def<python::object (const MoleculePerturber&, const RDKit::ROMol&, AtomIdx, std::mt19937&, const MolecularConstraints*)>(
      "InsertBond", InsertBond, (
      python::arg("molecule"),
      python::arg("begin_atom_idx"),
      python::arg("prng"),
      python::arg("constraints") = null_constraints))
    .def<python::object (const MoleculePerturber&, const RDKit::ROMol&, std::mt19937&, const MolecularConstraints*)>(
      "InsertBond", InsertBond, (
      python::arg("molecule"),
      python::arg("prng"),
      python::arg("constraints") = null_constraints))
    .def<void (MoleculePerturber::*)(MolecularPerturbationQueue&, const RDKit::ROMol&, AtomIdx, const MolecularConstraints*) const>(
      "BondInsertions", &MoleculePerturber::BondInsertions, (
      python::arg("queue"),
      python::arg("molecule"),
      python::arg("begin_atom_idx"),
      python::arg("constraints") = null_constraints))
    .def<void (MoleculePerturber::*)(MolecularPerturbationQueue&, const RDKit::ROMol&, const MolecularConstraints*) const>(
      "BondInsertions", &MoleculePerturber::BondInsertions, (
      python::arg("queue"),
      python::arg("molecule"),
      python::arg("constraints") = null_constraints))

    .def<python::object (const MoleculePerturber&, const RDKit::ROMol&, BondIdx, std::mt19937&, const MolecularConstraints*)>(
      "DeleteBond", DeleteBond, (
      python::arg("molecule"),
      python::arg("bond_idx"),
      python::arg("prng"),
      python::arg("constraints") = null_constraints))
    .def<python::object (const MoleculePerturber&, const RDKit::ROMol&, std::mt19937&, const MolecularConstraints*)>(
      "DeleteBond", DeleteBond, (
      python::arg("molecule"),
      python::arg("prng"),
      python::arg("constraints") = null_constraints))
    .def<void (MoleculePerturber::*)(MolecularPerturbationQueue&, const RDKit::ROMol&, BondIdx, const MolecularConstraints*) const>(
      "BondDeletions", &MoleculePerturber::BondDeletions, (
      python::arg("queue"),
      python::arg("molecule"),
      python::arg("bond_idx"),
      python::arg("constraints") = null_constraints))
    .def<void (MoleculePerturber::*)(MolecularPerturbationQueue&, const RDKit::ROMol&, const MolecularConstraints*) const>(
      "BondDeletions", &MoleculePerturber::BondDeletions, (
      python::arg("queue"),
      python::arg("molecule"),
      python::arg("constraints") = null_constraints))

    .def("Perturbation", Perturbation, (
      python::arg("molecule"),
      python::arg("prng"),
      python::arg("constraints") = null_constraints))
    .def("Perturbations", Perturbations, (
      python::arg("queue"),
      python::arg("molecule"),
      python::arg("constraints") = null_constraints))

    .def("EnablePerturbationType", &MoleculePerturber::EnablePerturbationType, (
      python::arg("perturbation_type"),
      python::arg("weight") = 1.0))
    .def("DisablePerturbationType", &MoleculePerturber::DisablePerturbationType, (
      python::arg("perturbation_type")))

    .def("GetPerturbationTypes", GetPerturbationTypes)
    .def("GetAtomicNumbers", GetAtomicNumbers)
    .def("GetRingAtomicNumbers", GetRingAtomicNumbers)
    .def("GetFormalCharges", GetFormalCharges)
    .def("GetNExplicitHydrogens", GetNExplicitHydrogens)
    .def("GetBondTypes", GetBondTypes)
    .def("GetRingBondTypes", GetRingBondTypes)

    .def("SetPerturbationTypes", SetPerturbationTypes, (
      python::arg("perturbation_types"),
      python::arg("weights"),
      python::arg("reset") = true))
    .def("SetAtomicNumbers", SetAtomicNumbers, (
      python::arg("atomic_numbers"),
      python::arg("weights")))
    .def("SetRingAtomicNumbers", SetRingAtomicNumbers, (
      python::arg("atomic_numbers"),
      python::arg("weights")))
    .def("SetFormalCharges", SetFormalCharges, (
      python::arg("formal_charges"),
      python::arg("weights")))
    .def("SetNExplicitHydrogens", SetNExplicitHydrogens, (
      python::arg("n_explicit_hydrogens"),
      python::arg("weights")))
    .def("SetBondTypes", SetBondTypes, (
      python::arg("bond_types"),
      python::arg("weights")))
    .def("SetRingBondTypes", SetRingBondTypes, (
      python::arg("bond_types"),
      python::arg("weights")))

    .def_readwrite("default_atomic_number", &MoleculePerturber::default_atomic_number)
    .def_readwrite("default_formal_charge", &MoleculePerturber::default_formal_charge)
    .def_readwrite("default_n_explicit_hydrogens", &MoleculePerturber::default_n_explicit_hydrogens)
    .def_readwrite("default_bond_type", &MoleculePerturber::default_bond_type)

    .def_readwrite("allow_disconnections", &MoleculePerturber::allow_disconnections)
    .def_readwrite("assess_connectivity_with_sssr", &MoleculePerturber::assess_connectivity_with_sssr)
    .def_readwrite("assess_distances_with_distance_matrix", &MoleculePerturber::assess_distances_with_distance_matrix)

    .def_readwrite("cyclicity_based_atomic_numbers", &MoleculePerturber::cyclicity_based_atomic_numbers)
    .def_readwrite("cyclicity_based_bond_types", &MoleculePerturber::cyclicity_based_bond_types)

    .def_readwrite("atom_insertion_max_n_neighbors", &MoleculePerturber::atom_insertion_max_n_neighbors)
    .def_readwrite("atom_insertion_drop_an_atom", &MoleculePerturber::atom_insertion_drop_an_atom)
    .def_readwrite("atom_insertion_dropping_is_optional", &MoleculePerturber::atom_insertion_dropping_is_optional)
    .def_readwrite("atom_insertion_only_pericentral_atoms_as_neighbor_candidates", &MoleculePerturber::atom_insertion_only_pericentral_atoms_as_neighbor_candidates)
    .def_readwrite("atom_insertion_min_distance_neighbor", &MoleculePerturber::atom_insertion_min_distance_neighbor)
    .def_readwrite("atom_insertion_max_distance_neighbor", &MoleculePerturber::atom_insertion_max_distance_neighbor)
    .def_readwrite("atom_insertion_randomize_atomic_number", &MoleculePerturber::atom_insertion_randomize_atomic_number)
    .def_readwrite("atom_insertion_randomize_formal_charge", &MoleculePerturber::atom_insertion_randomize_formal_charge)
    .def_readwrite("atom_insertion_randomize_n_explicit_hydrogens", &MoleculePerturber::atom_insertion_randomize_n_explicit_hydrogens)
    .def_readwrite("atom_insertion_iterate_atomic_numbers", &MoleculePerturber::atom_insertion_iterate_atomic_numbers)
    .def_readwrite("atom_insertion_iterate_formal_charges", &MoleculePerturber::atom_insertion_iterate_formal_charges)
    .def_readwrite("atom_insertion_iterate_n_explicit_hydrogens", &MoleculePerturber::atom_insertion_iterate_n_explicit_hydrogens)
    .def_readwrite("atom_insertion_iterate_bond_types", &MoleculePerturber::atom_insertion_iterate_bond_types)

    .def_readwrite("atom_deletion_allow_reconnections", &MoleculePerturber::atom_deletion_allow_reconnections)
    .def_readwrite("atom_deletion_only_consider_neighbors_for_reconnection", &MoleculePerturber::atom_deletion_only_consider_neighbors_for_reconnection)
    .def_readwrite("atom_deletion_preserve_bond_types_during_reconnection", &MoleculePerturber::atom_deletion_preserve_bond_types_during_reconnection)
    .def_readwrite("atom_deletion_min_distance_reconnection_atom", &MoleculePerturber::atom_deletion_min_distance_reconnection_atom)
    .def_readwrite("atom_deletion_max_distance_reconnection_atom", &MoleculePerturber::atom_deletion_max_distance_reconnection_atom)

    .def_readwrite("bond_insertion_min_distance_partner", &MoleculePerturber::bond_insertion_min_distance_partner)
    .def_readwrite("bond_insertion_max_distance_partner", &MoleculePerturber::bond_insertion_max_distance_partner)
    .def_readwrite("bond_insertion_max_atom_n_rings_membership", &MoleculePerturber::bond_insertion_max_atom_n_rings_membership)

    .def_readwrite("bond_insertion_randomize_bond_type", &MoleculePerturber::bond_insertion_randomize_bond_type)
    .def_readwrite("bond_insertion_iterate_bond_types", &MoleculePerturber::bond_insertion_iterate_bond_types)

    .def_readwrite("bond_deletion_allow_reroutes", &MoleculePerturber::bond_deletion_allow_reroutes)
    .def_readwrite("bond_deletion_preserve_bond_types_during_reroute", &MoleculePerturber::bond_deletion_preserve_bond_types_during_reroute)
    .def_readwrite("bond_deletion_min_distance_reroute", &MoleculePerturber::bond_deletion_min_distance_reroute)
    .def_readwrite("bond_deletion_max_distance_reroute", &MoleculePerturber::bond_deletion_max_distance_reroute);
};

#endif // !_PY_MOLECULE_PERTURBER_HPP_
