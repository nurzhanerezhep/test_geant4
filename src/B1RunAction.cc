//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
// #include "B1Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction()
: G4UserRunAction(),
  fEdep_1(0.),
  fEdep2_1(0.),
  fEdep_2(0.),
  fEdep2_2(0.)
{ 
  // add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray); 

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager1 = G4AccumulableManager::Instance();
  accumulableManager1->RegisterAccumulable(fEdep_1);
  accumulableManager1->RegisterAccumulable(fEdep2_1);

  G4AccumulableManager* accumulableManager2 = G4AccumulableManager::Instance();
  accumulableManager2->RegisterAccumulable(fEdep_2);
  accumulableManager2->RegisterAccumulable(fEdep2_2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run*)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager1 = G4AccumulableManager::Instance();
  accumulableManager1->Reset();
  G4AccumulableManager* accumulableManager2 = G4AccumulableManager::Instance();
  accumulableManager2->Reset();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables 
  G4AccumulableManager* accumulableManager1 = G4AccumulableManager::Instance();
  accumulableManager1->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep_1  = fEdep_1.GetValue();
  G4double edep2_1 = fEdep2_1.GetValue();
  
  G4double rms1 = edep2_1 - edep_1*edep_1/nofEvents;
  if (rms1 > 0.) rms1 = std::sqrt(rms1); else rms1 = 0.;

  const B1DetectorConstruction* detectorConstruction1
   = static_cast<const B1DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass1 = detectorConstruction1->GetScoringVolume1()->GetMass();
  G4double dose1 = edep_1/mass1;
  G4double rmsDose1 = rms1/mass1;

  G4double edep_2  = fEdep_2.GetValue();
  G4double edep2_2 = fEdep2_2.GetValue();

  G4double rms2 = edep2_2 - edep_2*edep_2/nofEvents;
  if (rms2 > 0.) rms2 = std::sqrt(rms2); else rms2 = 0.;

  const B1DetectorConstruction* detectorConstruction2
   = static_cast<const B1DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass2 = detectorConstruction2->GetScoringVolume2()->GetMass();
  G4double dose2= edep_2/mass2;
  G4double rmsDose2 = rms2/mass2;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const B1PrimaryGeneratorAction* generatorAction
   = static_cast<const B1PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }
        
  // Print
  //  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }
  
  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " Cumulated dose per run, in scoring volume 1: "
     << G4BestUnit(dose1,"Dose") << " rms = " << G4BestUnit(rmsDose1,"Dose")
     << G4endl
     << " Cumulated dose per run, in scoring volume 2: "
     << G4BestUnit(dose2,"Dose") << " rms = " << G4BestUnit(rmsDose2,"Dose")
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::AddEdep(G4double edep_1, G4double edep_2)
{
  fEdep_1  += edep_1;
  fEdep2_1 += edep_1*edep_1;
  fEdep_2  += edep_2;
  fEdep2_2 += edep_2*edep_2;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

