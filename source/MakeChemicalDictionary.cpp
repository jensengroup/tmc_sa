#include "ChemicalDictionary.hpp"
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/MolOps.h>
#include <set> // For storing skipped iterations

std::set<int> skipIterations = {
    7228, 7227, 7229, 7230, 7231,
    7232, 7233, 7234, 7235, 7236}; // Define iterations to skip

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

  std::cout << "Creating ChemicalDictionary with circular atomic environments "
            << "of radius " << environment_radius << std::endl;

  RDKit::ROMOL_SPTR molecule;
  int iteration = 1; // Counter for iterations
  while (!supplier.atEnd()) {
    std::cout << "At iteration" << iteration << std::endl;
    // Check if the current iteration should be skipped
    // if (skipIterations.find(iteration) != skipIterations.end()) {
    //   std::cout << "Skipping iteration " << iteration << std::endl;
    //   iteration++;
    //   continue; // Skip this iteration
    // }

    // std::cout << "Dictionary size at iteration " << iteration << ": "
    //       << dictionary.size() << std::endl;
    if (iteration == 2) {
      break;
    }

    // if (iteration % 4000 == 0) {
    //     std::ostringstream filename;
    //     filename << output_file_path << "_iter_" << iteration << ".dict";
    //     dictionary.Save(filename.str());
    //     ChemicalDictionary dictionary(environment_radius);
    //   // dictionary.Save(output_file_path);
    // };

    // supplier.next();
    // RDKit::ROMOL_SPTR molecule;
    // delete &molecule;
    molecule.reset(supplier.next());
    // std::cout << "Just reset" << iteration << std::endl;
    if (!molecule) {
      std::cout << "No mol" << iteration << std::endl;
      continue;
    };
    if (kekulize) {
      RDKit::RWMol kekulized_molecule(*molecule);
      if (!RDKit::MolOps::KekulizeIfPossible(kekulized_molecule)) {
        std::cout << "could not keku" << iteration << std::endl;
        continue;
      };
      dictionary.AddMolecule(kekulized_molecule);
    } else {
      dictionary.AddMolecule(*molecule);
    };

    iteration++;
  };
  dictionary.BuildPartialKeyDictionaries();

  dictionary.Save(output_file_path);

  return 0;
};
