// inspect_tree.cpp
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TBranch.h"
#include "TLeaf.h"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: inspect_tree <filelist.txt> [treeName]\n";
    return 1;
  }
  std::ifstream fin(argv[1]);
  if (!fin) { std::cerr << "Cannot open filelist\n"; return 1; }
  std::string first;
  for (std::string line; std::getline(fin, line); ) { if (!line.empty()) { first = line; break; } }
  if (first.empty()) { std::cerr << "Empty filelist\n"; return 1; }

  std::unique_ptr<TFile> f(TFile::Open(first.c_str()));
  if (!f || f->IsZombie()) { std::cerr << "Cannot open " << first << "\n"; return 1; }

  // List TTrees
  std::cout << "== Keys in file: " << first << "\n";
  f->GetListOfKeys()->Print();

  // Try provided tree or best guess
  std::string treeName = (argc>=3)? argv[2] : "events";
  TTree* t = dynamic_cast<TTree*>(f->Get(treeName.c_str()));
  if (!t) t = dynamic_cast<TTree*>(f->Get("Delphes"));
  if (!t) {
    std::cerr << "Could not find tree 'events' or 'Delphes'.\n";
    return 1;
  }
  std::cout << "\n== Using tree: " << t->GetName() << " with " << t->GetEntries() << " entries\n\n";

  // Print branches with type info
  TObjArray* branches = t->GetListOfBranches();
  for (int i=0;i<branches->GetEntries();++i) {
    auto br = (TBranch*)branches->At(i);
    std::cout << br->GetName() << "  |  " << br->GetTitle();
    // Show leaves and their types
    auto leaves = br->GetListOfLeaves();
    if (leaves) {
      for (int j=0;j<leaves->GetEntries();++j) {
        auto lf = (TLeaf*)leaves->At(j);
        std::cout << "\n    - leaf: " << lf->GetName()
                  << "  type: " << lf->GetTypeName()
                  << "  leaf count: " << (lf->GetLeafCount()? lf->GetLeafCount()->GetName(): "none");
      }
    }
    std::cout << "\n";
  }
  return 0;
}

