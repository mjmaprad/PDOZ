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
/// \file optical/OpNovice2/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
//#include "DetectorMessenger.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction()
{
 fExpHall_x = fExpHall_y = fExpHall_z = 50*cm;
  fScintillator_x    = 10.0*mm; 
  fScintillator_y    = 10.0*mm;
  fScintillator_z    = 10.0*mm;

  fScintillator = nullptr;

  fScintillatorMPT    = new G4MaterialPropertiesTable();
  fWorldMPT   = new G4MaterialPropertiesTable();
  fSurfaceMPT = new G4MaterialPropertiesTable();
  f2SurfaceMPT = new G4MaterialPropertiesTable();
//G4NistManager* man = G4NistManager::Instance();

////////////////////////////////////////////////////////////
//////


  fScintillator_LV  = nullptr;
  fWorld_LV = nullptr;

  fScintillatorMaterial  = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
  //fWorldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
 fWorldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  ////delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4Element* H = new G4Element("H", "H", 1., 1.01*g/mole);
  G4Element* C = new G4Element("C", "C", 6., 12.01*g/mole);
G4Material*    fGreaseMaterial= new G4Material("Silicone", 1.032*g/cm3, 2);
    fGreaseMaterial->AddElement(C,9);
    fGreaseMaterial->AddElement(H,10);


  //
  // ------------ Generate & Add Material Properties Table ------------
  //

  G4double photonEnergy[] =
  {2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
   2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
   2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
   2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
   2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
   2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
   2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
   3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
   3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
   3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);

  //--------------------------------------------------
  // Air
  //--------------------------------------------------

  G4double refractiveIndex[] =
  { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00};

  assert(sizeof(refractiveIndex) == sizeof(photonEnergy));

  G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX", photonEnergy, refractiveIndex, nEntries);




  G4double absClad[] =
  {20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m};

  assert(sizeof(absClad) == sizeof(photonEnergy));

  //--------------------------------------------------
  // Silicone
  //--------------------------------------------------

   G4double refractiveIndexSilicone[] =
   { 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46,
     1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46,
     1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46,
     1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46,
     1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46};

   assert(sizeof(refractiveIndexSilicone) == sizeof(photonEnergy));
fSurface_scint_grease = new G4OpticalSurface("Surface_scint_grease");
fSurface_scint_grease ->SetType(dielectric_dielectric);
fSurface_scint_grease ->SetModel(unified);
fSurface_scint_grease -> SetSigmaAlpha(0.1);
fSurface_scint_grease ->SetFinish(ground);
  // Add entries into properties table
  G4MaterialPropertiesTable* mptSilicone = new G4MaterialPropertiesTable();
  mptSilicone->
           AddProperty("RINDEX",photonEnergy,refractiveIndexSilicone,nEntries);
  mptSilicone->AddProperty("ABSLENGTH",photonEnergy,absClad,nEntries);
fSurface_scint_grease->SetMaterialPropertiesTable(mptSilicone);
  fGreaseMaterial->SetMaterialPropertiesTable(mptSilicone);




/////////////////////////////////////////////////////////////////////////////////////////////////////7/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4NistManager* man = G4NistManager::Instance();

G4Material *fSiPMMaterial = man->FindOrBuildMaterial("G4_Si");
G4MaterialPropertiesTable* silicon = new G4MaterialPropertiesTable();

fSurface_grease_sipm = new G4OpticalSurface("fSurface__grease_sipm");
fSurface_grease_sipm ->SetType(dielectric_metal);
fSurface_grease_sipm ->SetModel(unified);
fSurface_grease_sipm ->SetFinish(ground);
fSurface_grease_sipm -> SetSigmaAlpha(0.0);
fSurface_grease_sipm->SetMaterialPropertiesTable(silicon);
fSiPMMaterial->SetMaterialPropertiesTable(silicon);



//G4Material *fPb = man->FindOrBuildMaterial("G4_Pb");

//G4Material *fTeflon = man->FindOrBuildMaterial("G4_MYLAR");




 //TiO2
 ////https://refractiveindex.info/?shelf=main&book=TiO2&page=Siefke
//sigma_alpha = 0.1;
fSurface_tio2 = new G4OpticalSurface("fSurface_tio2");
fSurface_tio2->SetType(dielectric_dielectric);
fSurface_tio2->SetModel(unified);
fSurface_tio2->SetFinish(groundbackpainted);
//fSurface_tio2 -> SetSigmaAlpha(sigma_alpha);
G4double Energy_tio2[5]={0.4*eV,2.5*eV,5*eV,7.5*eV,10*eV};
G4double refin_tio2[5]={2.22,2.5,2.37,1.73,1.11};
//S. Yin, M. Komatsu, Q. Zhang, F. Saito, T. Sato, J. Mater. Sci., 42(7), 2399 (2007). 
//SIZE DEPENDENT REFLECTIVE PROPERTIES OF TIO2 NANOPARTICLES AND REFLECTORS MADE THEREOF 
G4double reflect_tio2[5]={0.91,0.91,0.91,0.91,0.91};
 f2SurfaceMPT->AddProperty("RINDEX",Energy_tio2,refin_tio2,2);
 f2SurfaceMPT->AddProperty("REFLECTIVITY",Energy_tio2,reflect_tio2,2);
fSurface_tio2->SetMaterialPropertiesTable(f2SurfaceMPT);
G4Material *ftio2 = man->FindOrBuildMaterial("G4_TITANIUM_DIOXIDE");
ftio2->SetMaterialPropertiesTable(f2SurfaceMPT);

////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //--------------------------------------------------
  // CSI_TL
  //--------------------------------------------------
  G4Element* eCs= new G4Element("Caesium", "Cs", 55, 132.9055*g/mole);
  G4Element* eI = new G4Element("Iodine", "I", 53, 126.9044*g/mole);
  G4Material *fCSI= new G4Material("fCSI", 4.51*g/cm3, 2);
  fCSI->AddElement(eCs,1);
  fCSI->AddElement(eI,1);
G4Material *fTl=G4NistManager::Instance()->FindOrBuildMaterial("G4_Tl");
////https://indico.cern.ch/event/771235/contributions/3205236/attachments/1778654/2894242/geant4_spdak2019_hanwook_2.pdf
  G4Material *fCSI_TL= new G4Material("fCSI_TL", 4.51*g/cm3, 2);
  fCSI_TL->AddMaterial(fTl,0.001*perCent);
  fCSI_TL->AddMaterial(fCSI,(100-0.001)*perCent);

//https://www.crystals.saint-gobain.com/sites/imdf.crystals.com/files/documents/csitl-and-na-material-data-sheet.pdf
  G4double EnergyCSI_TL[] =
    {3.11*eV,2.80*eV,2.55*eV,2.41*eV,2.27*eV,2.11*eV,1.97*eV,1.87*eV,1.75*eV

    };

  const G4int nEntriesCSI_TL = sizeof(EnergyCSI_TL)/sizeof(G4double);

  G4double refractiveIndexCSI_TL[] =
    {1.79,     1.79,     1.79,     1.79,     1.79,     1.79,     1.79,     1.79,     1.79
    };

  assert(sizeof(refractiveIndexCSI_TL) == sizeof(EnergyCSI_TL));


  G4double absCSI_TL[] =
    {10*cm,     10*cm,     10*cm,	
     10*cm,     10*cm,     10*cm,	
     10*cm,     10*cm,     10*cm
    };

  assert(sizeof(absCSI_TL) == sizeof(EnergyCSI_TL));

  G4double scintillationFastCSI_TL[] =
    {
 1.13E-01,2.42E-01,5.91E-01,8.49E-01,9.84E-01,7.84E-01,5.26E-01,3.01E-01,1.24E-01

    }; 

  assert(sizeof(scintillationFastCSI_TL) == sizeof(EnergyCSI_TL));
  /// In CSI_TL the emission wavelength peak is 520 nm

  // Add entries into properties table
  G4MaterialPropertiesTable* mptCSI_TL = new G4MaterialPropertiesTable();
  mptCSI_TL->AddProperty("RINDEX",EnergyCSI_TL,refractiveIndexCSI_TL,nEntriesCSI_TL);
  mptCSI_TL->AddProperty("ABSLENGTH",EnergyCSI_TL,absCSI_TL,nEntriesCSI_TL);
  mptCSI_TL->AddProperty("FASTCOMPONENT",EnergyCSI_TL, scintillationFastCSI_TL,nEntriesCSI_TL);
  //mptCSI_TL->AddProperty("SLOWCOMPONENT",EnergyCSI_TL, scintillationSlowCSI_TL,nEntriesCSI_TL);
  mptCSI_TL->AddConstProperty("FASTTIMECONSTANT", 1000*ns);///
  mptCSI_TL->AddConstProperty("SLOWTIMECONSTANT", 3500*ns);//
  mptCSI_TL->AddConstProperty("SCINTILLATIONYIELD",54000./MeV);
  mptCSI_TL->AddConstProperty("YIELDRATIO",1.0);
  mptCSI_TL->AddConstProperty("RESOLUTIONSCALE",1);//
fCSI_TL->SetMaterialPropertiesTable(mptCSI_TL);
// Set the Birks Constant for the BC408 scintillator
fCSI_TL->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////////////////////7/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //--------------------------------------------------
  //  BC408
  // //--------------------------------------------------
  // G4Element* H = new G4Element("H", "H", 1., 1.01*g/mole);
  // G4Element* C = new G4Element("C", "C", 6., 12.01*g/mole);
    fBC408= new G4Material("BC408", 1.032*g/cm3, 2);
    fBC408->AddElement(C,9);
    fBC408->AddElement(H,10);

G4double EnergyBC408[] = 
{3.299*eV, 3.201*eV, 3.130*eV, 3.089*eV, 3.069*eV, 3.056*eV, 3.036*eV, 3.023*eV, 3.016*eV, 3.004*eV, 2.991*eV, 2.972*eV, 2.947*eV, 2.923*eV, 2.905*eV, 2.887*eV, 2.863*eV, 2.852*eV, 2.829*eV, 2.818*eV, 2.818*eV, 2.795*eV, 2.768*eV, 2.736*eV, 2.699*eV, 2.658*eV, 2.638*eV, 2.599*eV, 2.562*eV, 2.511*eV, 2.455*eV, 2.421*eV, 2.384*eV};
const G4int nEntriesBC408 = sizeof(EnergyBC408)/sizeof(G4double);


G4double refractiveIndexBC408[] =
{ 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58};
assert(sizeof(refractiveIndexBC408) == sizeof(EnergyBC408));


G4double absBC408[] =
{210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm};
assert(sizeof(absBC408) == sizeof(EnergyBC408));


G4double scintillationFastBC408[] =
{0.057, 0.097, 0.163, 0.260, 0.326, 0.396, 0.449, 0.520, 0.564, 0.648, 0.714, 0.793, 0.877, 0.943, 0.978, 0.996, 0.965, 0.930, 0.833, 0.753, 0.705, 0.599, 0.537, 0.476, 0.388, 0.282, 0.242, 0.181, 0.141, 0.088, 0.048, 0.035, 0.022}; 
assert(sizeof(scintillationFastBC408) == sizeof(EnergyBC408));

// Add entries into properties table
G4MaterialPropertiesTable* BC408 = new G4MaterialPropertiesTable();
BC408->AddProperty("RINDEX",EnergyBC408,refractiveIndexBC408,nEntriesBC408);
BC408->AddProperty("ABSLENGTH",EnergyBC408,absBC408,nEntriesBC408);
BC408->AddProperty("FASTCOMPONENT",EnergyBC408, scintillationFastBC408,nEntriesBC408);
//BC408->AddProperty("SLOWCOMPONENT",EnergyBC408, scintillationSlowBC408,nEntriesBC408);
BC408->AddConstProperty("FASTTIMECONSTANT", 0.9*ns);
BC408->AddConstProperty("SLOWTIMECONSTANT", 2.1*ns);
BC408->AddConstProperty("SCINTILLATIONYIELD",10240./MeV);
BC408->AddConstProperty("YIELDRATIO",1.0);
BC408->AddConstProperty("RESOLUTIONSCALE",1.0);
fBC408->SetMaterialPropertiesTable(BC408);
// Set the Birks Constant for the BC408 scintillator
fBC408->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // ---------------------- Volumes ------------------------
  // The experimental Hall
  G4Box* world_box = new G4Box("Worlds", fExpHall_x, fExpHall_y, fExpHall_z);

  fWorld_LV
    = new G4LogicalVolume(world_box, fWorldMaterial, "fWorld_LV", 0, 0, 0);
G4cout<<"the world is "<<fWorldMaterial->GetName()<<G4endl;

  G4VPhysicalVolume* world_PV
    = new G4PVPlacement(0, G4ThreeVector(), fWorld_LV, "world_PV", 0, false, 0,true);
//-----------------------------------------------------------
/////////////////////////////////////////////////
G4double Opt_win=2*mm;
  	G4double innerRadiusOfTheTube = 0. * mm;///// 
	G4double outerRadiusOfTheTube = 3 * mm;// x and y of grease
	G4double hightOfTheTube = 0.5 * mm;///height of sipm
	G4double hightOfTheTube2 = 5 * um;//////z of grease
	G4double startAngleOfTheTube = 0. * deg;
	G4double spanningAngleOfTheTube = 360. * deg;
	 G4double gapBetweenHoles = (fScintillator_x-outerRadiusOfTheTube*2)/3+outerRadiusOfTheTube;
  ///////////////////////////////////////////////
  //teflon
  G4double tef_thick=100*um;
 G4double mytefX = fScintillator_x+tef_thick*2;
 G4double mytefY =fScintillator_y+tef_thick*2;
G4double mytefZ =fScintillator_z+tef_thick; 
// G4double	mytefz = fScintillator_z5.1 * mm;
	///////////////////// Grease Geometry//////////////////////////////
	// G4Tubs* SolidGrease = new G4Tubs("Grease", innerRadiusOfTheTube,
	// 		outerRadiusOfTheTube, hightOfTheTube2 / 2, startAngleOfTheTube,
	// 		spanningAngleOfTheTube);
  G4Box* SolidOpt_win = new G4Box("SolidOpt_win",
			fScintillator_x / 2, fScintillator_y / 2,Opt_win / 2);
	fLogiOpt_win = new G4LogicalVolume(SolidOpt_win, fGreaseMaterial, "fLogiOpt_win");
/////////////////////////////////////////////////////////////////////////////
	// fLogicGrease = new G4LogicalVolume(SolidGrease, fGreaseMaterial, "Grease");
  G4Box* SolidGrease = new G4Box("Grease",
			outerRadiusOfTheTube / 2, outerRadiusOfTheTube / 2,hightOfTheTube2 / 2);
	fLogicGrease = new G4LogicalVolume(SolidGrease, fGreaseMaterial, "Grease");
// 	///////////////////// end - Grease Geometry//////////////////////////////

// 	///////////////////// SiPM Geometry//////////////////////////////
	// G4Tubs* SolidSiPM = new G4Tubs("SiPM", innerRadiusOfTheTube, outerRadiusOfTheTube,
	// 			hightOfTheTube / 2, startAngleOfTheTube, spanningAngleOfTheTube);

	// 	fLogicSiPM = new G4LogicalVolume(SolidSiPM, fSiPMMaterial, "SiPM");
  G4Box* SolidSiPM = new G4Box("SiPM",
			outerRadiusOfTheTube / 2, outerRadiusOfTheTube / 2,hightOfTheTube / 2);




		fLogicSiPM = new G4LogicalVolume(SolidSiPM, fSiPMMaterial, "SiPM");
    /////////////////////////////////////////////////////////////////////////
    G4double Dis_tra_scin=2.5*mm+fScintillator_y/2; ///this 2.5 is from the center so the distance is 5mm
// 	///////////////////// Grease - SiPM Geometry//////////////////////////////
 	G4int nOfSolids = 2;
 	for (G4int k = 1; k <= nOfSolids; k++) {

		G4String volNameGrease;
		G4String volNameSiPM;
		switch (k) {
		case 1:
			//volNameGrease = "FOpt_winPV";
			volNameSiPM = "SiPM1";
      //
 		FOpt_winPV = new G4PVPlacement(0, G4ThreeVector(0,Dis_tra_scin,-fScintillator_z/2-Opt_win/2),fLogiOpt_win,"FOpt_winPV", fWorld_LV, false, 0,true);
    fsipmPV=new G4PVPlacement(0, G4ThreeVector(-gapBetweenHoles/2,Dis_tra_scin,-fScintillator_z/2-Opt_win-hightOfTheTube/2),fLogicSiPM,volNameSiPM, fWorld_LV, false, 0,true);
			break;
		case 2:
			//volNameGrease = "Grease2";
			volNameSiPM = "SiPM2";
     //
    // fgreasePV = new G4PVPlacement(0, G4ThreeVector(+gapBetweenHoles/2,Dis_tra_scin,-fScintillator_z/2-hightOfTheTube2/2),fLogicGrease,volNameGrease, fWorld_LV, false, 0,true);
    fsipmPV=new G4PVPlacement(0, G4ThreeVector(+gapBetweenHoles/2,Dis_tra_scin,-fScintillator_z/2-Opt_win-hightOfTheTube/2),fLogicSiPM,volNameSiPM, fWorld_LV, false, 0,true);

			break;
			}
   }
   ////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////
   /////////////////////////////////



G4Box* solidScintillator = new G4Box("solidScintillator",
			fScintillator_x / 2, fScintillator_y / 2,fScintillator_z / 2);
	fScintillator_LV = new G4LogicalVolume(solidScintillator,fCSI_TL,
			"fScintillator_LV");
G4cout<<"the big scintillator is "<<fCSI_TL->GetName()<<G4endl;

	fScintillator = new G4PVPlacement(0, G4ThreeVector(0., Dis_tra_scin, 0.),
			fScintillator_LV, "Scintillator", fWorld_LV, false, 0,true);
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/////////////////////////// tio2
  G4double tio_thick=300*um;
 G4double mytioX = fScintillator_x+tio_thick*2;
 G4double mytioY =fScintillator_y+tio_thick*2;
G4double mytioZ =fScintillator_z+tio_thick+Opt_win; 
G4Box* solidtio = new G4Box("solidtio",
			mytioX / 2, mytioY / 2,mytioZ / 2);
      G4Box* solidtioHH = new G4Box("solidtioHH",
			fScintillator_x / 2, fScintillator_y / 2,(fScintillator_z+Opt_win) / 2);
	G4SubtractionSolid* subtracttio1 = new G4SubtractionSolid("subtracttio1",
			solidtio, solidtioHH, 0,
			G4ThreeVector(
					0,
					0,-tio_thick/2));


	// G4SubtractionSolid* subtracttio2 = new G4SubtractionSolid("subtracttio2",
	// 		subtracttio1, SolidOpt_win, 0,
	// 		G4ThreeVector(
	// 				0,
	// 				0,-tio_thick/2-Opt_win)); 


	G4LogicalVolume* logictio=	 new G4LogicalVolume(subtracttio1, //its solid
						             ftio2,    //its material
							     "logictio");           //its name

	 G4PVPlacement* ftio = new G4PVPlacement(0,                             //no rotation
				G4ThreeVector(0,Dis_tra_scin,tio_thick/2-Opt_win/2),                  //at (0,0,0)
				logictio,                      //its logical volume
				"ftio",                           //its name
				fWorld_LV,                                //its mother  volume
				false,                            //no boolean operation
				0,true);                               //copy number
/////////////////////////// Al cover
  G4double al_thick=1.5*mm;
   G4double myalX = mytioX+al_thick*2;
 G4double myalY =mytioY+al_thick*2;
G4double myalZ =mytioZ+al_thick; 
G4Material *falum = man->FindOrBuildMaterial("G4_Al");

G4Box* solidal = new G4Box("solidal",
			myalX / 2, myalY / 2,myalZ / 2);
	G4SubtractionSolid* subtractal1 = new G4SubtractionSolid("subtractal1",
			solidal, solidtio, 0,
			G4ThreeVector(
					0,
					0,-al_thick/2));
	G4LogicalVolume* logical=	 new G4LogicalVolume(subtractal1, //its solid
						             falum,    //its material
							     "logical");           //its name

	G4PVPlacement * fal = new G4PVPlacement(0,                             //no rotation
				G4ThreeVector(0,Dis_tra_scin,+al_thick/2+tio_thick/2-Opt_win/2),                  //at (0,0,0)
				logical,                      //its logical volume
				"fal",                           //its name
				fWorld_LV,                                //its mother  volume
				false,                            //no boolean operation
				0,true);                               //copy number

/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//// the lower dosimeter
///////////////////////////////////////////////////////////////////////////////////


 //Teflon-scint
G4double sigma_alpha = 0.1;
fSurface = new G4OpticalSurface("fSurface");
fSurface->SetType(dielectric_dielectric);
fSurface->SetModel(unified);
fSurface->SetFinish(groundbackpainted);
fSurface -> SetSigmaAlpha(sigma_alpha);
G4double Energy_teflon[2]={2.755*eV,2.755*eV};
G4double refin_teflon[2]={1.35,1.35};
//G4double transmit_teflon[2]={0.06,0.06};
G4double reflect_teflon[2]={0.95,0.95};
 fSurfaceMPT->AddProperty("RINDEX",Energy_teflon,refin_teflon,2);
 fSurfaceMPT->AddProperty("REFLECTIVITY",Energy_teflon,reflect_teflon,2);
fSurface->SetMaterialPropertiesTable(fSurfaceMPT);
G4Material *fTeflon = man->FindOrBuildMaterial("G4_TEFLON");
fTeflon->SetMaterialPropertiesTable(fSurfaceMPT);
///////////
////////////////////////////////////////////////////////////////

 G4double fScintillator_L_z= 1.5*mm;

for (G4int k = 1; k <= nOfSolids; k++) {

		G4String volNameGrease;
		G4String volNameSiPM;
		switch (k) {
		case 1:
			volNameGrease = "Grease3";
			volNameSiPM = "SiPM3";
      //
 		fgreasePV = new G4PVPlacement(0, G4ThreeVector(-gapBetweenHoles/2,-Dis_tra_scin,+ myalZ/2-Opt_win/2-fScintillator_L_z-hightOfTheTube2/2+(al_thick/2+tio_thick/2-tef_thick)),fLogicGrease,volNameGrease, fWorld_LV, false, 0,true);
    fsipmPV=new G4PVPlacement(0, G4ThreeVector(-gapBetweenHoles/2,-Dis_tra_scin,+ myalZ/2-Opt_win/2-fScintillator_L_z-hightOfTheTube2-hightOfTheTube/2+(al_thick/2+tio_thick/2-tef_thick)),fLogicSiPM,volNameSiPM, fWorld_LV, false, 0,true);
			break;
		case 2:
			volNameGrease = "Grease4";
			volNameSiPM = "SiPM4";
     //
     fgreasePV = new G4PVPlacement(0, G4ThreeVector(+gapBetweenHoles/2,-Dis_tra_scin,+ myalZ/2-Opt_win/2-fScintillator_L_z-hightOfTheTube2/2+(al_thick/2+tio_thick/2-tef_thick)),fLogicGrease,volNameGrease, fWorld_LV, false, 0,true);
    fsipmPV=new G4PVPlacement(0, G4ThreeVector(+gapBetweenHoles/2,-Dis_tra_scin,+ myalZ/2-Opt_win/2-fScintillator_L_z-hightOfTheTube2-hightOfTheTube/2+(al_thick/2+tio_thick/2-tef_thick)),fLogicSiPM,volNameSiPM, fWorld_LV, false, 0,true);

			break;
			}
   }

   //////
   G4Box* solidScintillator_L = new G4Box("solidScintillator_L",
			fScintillator_x / 2, fScintillator_y / 2,fScintillator_L_z / 2);
	fScintillator_L_LV = new G4LogicalVolume(solidScintillator_L, fBC408,
			"fScintillator_L_LV");
G4cout<<"the small scintillator is "<<fBC408->GetName()<<G4endl;
fScintillator_L = new G4PVPlacement(0, G4ThreeVector(0., -Dis_tra_scin, myalZ/2-Opt_win/2-fScintillator_L_z/2+al_thick/2+tio_thick/2-tef_thick),
			fScintillator_L_LV, "Scintillator_PV", fWorld_LV, false, 0,true);

/////////////////////////
G4double mytefZ_L =fScintillator_L_z+tef_thick;
G4Box* solidTeflon_L = new G4Box("solidTeflon_L",
			mytefX / 2, mytefY / 2,mytefZ_L / 2);
G4SubtractionSolid* subtracttef2 = new G4SubtractionSolid("subtracttef2",
			solidTeflon_L, solidScintillator_L, 0,
			G4ThreeVector(
					0,
					0,-tef_thick/2));

	G4LogicalVolume* logicTeflon_2=	 new G4LogicalVolume(subtracttef2, //its solid
						             fTeflon,    //its material
							     "logicTeflon_2");           //its name

	G4VPhysicalVolume* fFoil_LV= new G4PVPlacement(0,                             //no rotation
				G4ThreeVector(0,-Dis_tra_scin,myalZ/2-Opt_win/2-fScintillator_L_z/2+al_thick/2+tio_thick/2-tef_thick/2),                  //at (0,0,0)
				logicTeflon_2,                      //its logical volume
				"fFoil_LV",                           //its name
				fWorld_LV,                                //its mother  volume
				false,                            //no boolean operation
				0,true);                               //copy number


/////////////////////////////////////
/////////////////////////////////////// box around the scints
/////////////////////////////////////
 G4Material *fplexyglass = man->FindOrBuildMaterial("G4_PLEXIGLASS");

  G4double color_thick=2*mm;
 G4double mycolorX = Dis_tra_scin*2;
 G4double mycolorY =Dis_tra_scin*4;
G4double mycolorZ =myalZ; 
G4cout<<"the PMMA dimensions: x = "<<mycolorX+ color_thick<<"y = "<<mycolorY/ 2+ color_thick<<"Z = "<<mycolorZ / 2+ color_thick/2<<G4endl;
G4Box* solidcolor = new G4Box("solidcolor",
			mycolorX/2+ color_thick,mycolorY/ 2+ color_thick,mycolorZ / 2+ color_thick/2);

      G4Box* solidcolor2 = new G4Box("solidcolor2",
			mycolorX/2,mycolorY/ 2,mycolorZ / 2);

      G4SubtractionSolid* subtractcolor1 = new G4SubtractionSolid("subtractcolor1",
			solidcolor, solidcolor2, 0,
			G4ThreeVector(
					0,
					0,-color_thick/2));
// 	
	G4LogicalVolume* logiccolor=	 new G4LogicalVolume(subtractcolor1, //its solid
						             fplexyglass,    //its material
							     "logiccolor");           //its name

	new G4PVPlacement(0,                             //no rotation
				G4ThreeVector(0,0,+color_thick/2+tio_thick/2+al_thick/2-Opt_win/2),                  //at (0,0,0)
				logiccolor,                      //its logical volume
				"color_pv",                           //its name
				fWorld_LV,                                //its mother  volume
				false,                            //no boolean operation
				0,true);                               //copy number
///////////////////////////////////////////////////////////////////////////////////////// outer box
G4double L_L=100*mm;
G4Box* solidoutbox = new G4Box("solidoutbox",
			mycolorX/2+ color_thick+L_L,mycolorY/ 2+ color_thick+L_L,mycolorZ / 2+ color_thick/2+L_L/2);
           
            G4SubtractionSolid* subtractoutbox1 = new G4SubtractionSolid("subtractoutbox1",
			solidoutbox, solidcolor, 0,
			G4ThreeVector(
					0,
					0,-L_L/2));

	G4LogicalVolume* logioutbox=	 new G4LogicalVolume(subtractoutbox1, //its solid
						             fWorldMaterial,    //its material
							     "logioutbox");           //its name

	new G4PVPlacement(0,                             //no rotation
				G4ThreeVector(0,0,L_L/2+color_thick/2+tio_thick/2+al_thick/2-Opt_win/2),                  //at (0,0,0)
				logioutbox,                      //its logical volume
				"outboxpv_pv",                           //its name
				fWorld_LV,                                //its mother  volume
				false,                            //no boolean operation
				0,true);                               //copy number
// /////////////////////////////////////////////////////////////////////////////////////////



  // ------------- Surface --------------
 
   G4LogicalBorderSurface* surface =
           new G4LogicalBorderSurface("Surface",
                                  fScintillator_L, fFoil_LV, fSurface);
 G4LogicalBorderSurface* surfacetio2 =
           new G4LogicalBorderSurface("surfacetio2",
                                  fScintillator, ftio, fSurface_tio2);

  //  G4LogicalBorderSurface* surface1 =
  //          new G4LogicalBorderSurface("Surface1",
  //                                 fScintillator, hole1_PV, fSurface_world);
  //    G4LogicalBorderSurface* surface2 =
  //          new G4LogicalBorderSurface("Surface2",
  //                                 fScintillator, hole2_PV, fSurface_world);   

  //  G4LogicalBorderSurface* surface3 =
  //          new G4LogicalBorderSurface("Surface3",
  //                                  hole1_PV,fgreasePV, fSurface_world);

// G4LogicalBorderSurface* surfacetio3 =
//            new G4LogicalBorderSurface("surfacetio3",
//                                   ftio, FOpt_winPV, fSurface_scint_grease);

      G4LogicalBorderSurface* surface7=
           new G4LogicalBorderSurface("Surface7",
                                   fScintillator,FOpt_winPV, fSurface_scint_grease);     
     G4LogicalBorderSurface* surface4 =
           new G4LogicalBorderSurface("Surface4",
                                   FOpt_winPV,fsipmPV, fSurface_grease_sipm);   
                                                                                              
      G4LogicalBorderSurface* surface5 =
           new G4LogicalBorderSurface("Surface5",
                                   fScintillator_L,fgreasePV, fSurface_scint_grease);    
      G4LogicalBorderSurface* surface6 =
           new G4LogicalBorderSurface("Surface6",
                                   fgreasePV,fsipmPV, fSurface_grease_sipm);    



  return world_PV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSurfaceSigmaAlpha(G4double v) {
  // fSurface->SetSigmaAlpha(v);
  // fSurface_scint_grease->SetSigmaAlpha(0);
  // fSurface_grease_sipm->SetSigmaAlpha(0);
  // G4RunManager::GetRunManager()->GeometryHasBeenModified();

  // G4cout << "Surface sigma alpha set to: " << fSurface->GetSigmaAlpha()
  //        << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSurfacePolish(G4double v) {
  fSurface->SetPolish(v);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4cout << "Surface polish set to: " << fSurface->GetPolish()
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddScintillatorMPV(const char* c,
                                     G4MaterialPropertyVector* mpv) {
  fScintillatorMPT->AddProperty(c, mpv);
  G4cout << "The MPT for the box is now: " << G4endl;
  fScintillatorMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddWorldMPV(const char* c,
                                       G4MaterialPropertyVector* mpv) {
  fWorldMPT->AddProperty(c, mpv);
  G4cout << "The MPT for the world is now: " << G4endl;
  fWorldMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddSurfaceMPV(const char* c,
                                         G4MaterialPropertyVector* mpv) {
  fSurfaceMPT->AddProperty(c, mpv);
  G4cout << "The MPT for the surface is now: " << G4endl;
  fSurfaceMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddScintillatorMPC(const char* c, G4double v) {
  fScintillatorMPT->AddConstProperty(c, v);
  G4cout << "The MPT for the box is now: " << G4endl;
  fScintillatorMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddWorldMPC(const char* c, G4double v) {
  fWorldMPT->AddConstProperty(c, v);
  G4cout << "The MPT for the world is now: " << G4endl;
  fWorldMPT->DumpTable();
  G4cout << "............." << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddSurfaceMPC(const char* c, G4double v) {
  fSurfaceMPT->AddConstProperty(c, v);
  G4cout << "The MPT for the surface is now: " << G4endl;
  fSurfaceMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetWorldMaterial(const G4String& mat) {
  G4Material* pmat = G4NistManager::Instance()->FindOrBuildMaterial(mat);
  if (pmat && fWorldMaterial != pmat) {
    fWorldMaterial = pmat;
    if (fWorld_LV) {
      fWorld_LV->SetMaterial(fWorldMaterial);
      fWorldMaterial->SetMaterialPropertiesTable(fWorldMPT);
    }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    G4cout << "World material set to " << fWorldMaterial->GetName()
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetScintillatorMaterial(const G4String& mat) {
  G4Material* pmat = G4NistManager::Instance()->FindOrBuildMaterial(mat);
  if (pmat && fScintillatorMaterial != pmat) {
    fScintillatorMaterial = pmat;
    if (fScintillator_LV) {
      fScintillator_LV->SetMaterial(fScintillatorMaterial);
      fScintillatorMaterial->SetMaterialPropertiesTable(fScintillatorMPT);
      fScintillatorMaterial->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
    }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    G4cout << "Scintillator material set to " << fScintillatorMaterial->GetName()
           << G4endl;
  }
}