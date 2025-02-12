#include "ChemicalDictionary.hpp"
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/MolOps.h>

int main(int argc, char *argv[]) {
  if (argc < 3) {
    return 1;
  };

  bool kekulize = false;
  unsigned environment_radius = 2;
  std::string input_file_path(argv[1]);
  std::string output_file_path(argv[2]);
  if (argc > 3) {
    environment_radius = std::stoi(argv[3]);
  };

  RDKit::SmilesMolSupplier supplier(input_file_path, " \t", 0, 1, false, true);

  ChemicalDictionary dictionary(environment_radius);

  std::cout
      << "Creating ChemicalDictionary object with circular atomic (radius "
      << environment_radius
      << "). This object contains only the first entry of the given .smi "
         "file.\n"
         "Any remaining TMC SMILES should be added using the Python "
         "functionality."
      << " Radius of environment: " << environment_radius << std::endl;

  RDKit::ROMOL_SPTR molecule;
  int iteration = 1; // Counter for iterations
  while (!supplier.atEnd()) {
    std::cout << "Loading molecule " << iteration << std::endl;
    molecule.reset(supplier.next());
    if (!molecule) {
      std::cout << "No mol obtained from the smiles at line " << iteration
                << " . Checking next molecule" << std::endl;
      continue;
    };
    if (kekulize) {
      RDKit::RWMol kekulized_molecule(*molecule);
      if (!RDKit::MolOps::KekulizeIfPossible(kekulized_molecule)) {
        std::cout << "Could not kekulize" << iteration << std::endl;
        continue;
      };
      dictionary.AddMolecule(kekulized_molecule);
    } else {
      dictionary.AddMolecule(*molecule);
    };

    // We break as soon as we have molecule to create a single entry dict at the
    // c++ level. The current c++ implementation can not handle the complexity
    // of TMC fragments memory wise. This single entry dict can be loaded in
    // Python where the memory issue is not present.
    if (molecule) {
      std::cout
          << "Breaking after 1 iteration. Add the remaining SMILES with Python"
          << std::endl;
      break;
    }
    iteration++;
  };
  dictionary.BuildPartialKeyDictionaries();

  dictionary.Save(output_file_path);

  return 0;
};
