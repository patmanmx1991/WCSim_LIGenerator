//  -*- mode:c++; tab-width:4;  -*-
#include "WCSimDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4IntersectionSolid.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4PVReplica.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4AssemblyVolume.hh"
#include "G4SubtractionSolid.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalSurface.hh"
#include "G4UserLimits.hh"
#include "G4ReflectionFactory.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "WCSimTuningParameters.hh" //jl145

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4VSolid.hh" // following 4 .hh required for CADMesh.
#include "G4Trd.hh"
#include "G4NistManager.hh"
#include "G4AssemblyVolume.hh"


struct RelasticPlacementConfiguration {

    G4ThreeVector EntireDetectOffset; //< Relative offset for

    G4Material* InnerDetectorMaterial;
    G4Material* OuterDetectorMaterial;

    G4VisAttributes* BlackTyvekVis;
    G4Material* BlackTyvekMaterial;
    G4double BlackTyvekInnerRadius;
    G4double BlackTyvekOuterRadius;
    G4double BlackTyvekBarrelLength;

    G4VisAttributes* DeadSpaceVis;
    G4Material* DeadSpaceMaterial;
    G4double DeadSpaceInnerRadius;
    G4double DeadSpaceOuterRadius;
    G4double DeadSpaceBarrelLength;

    G4VisAttributes* OuterDetectorVis;
    G4double OuterDetectorInnerRadius;
    G4double OuterDetectorOuterRadius;
    G4double OuterDetectorBarrelLength;

    G4VisAttributes* InnerDetectorVis;
    G4double InnerDetectorInnerRadius;
    G4double InnerDetectorOuterRadius;
    G4double InnerDetectorBarrelLength;

    G4VisAttributes* WhiteTyvekVis;
    G4Material* WhiteTyvekMaterial;
    G4double WhiteTyvekInnerRadius;
    G4double WhiteTyvekOuterRadius;
    G4double WhiteTyvekBarrelLength;

    G4VisAttributes* WallTyvekVis;
    G4Material* WallTyvekMaterial;
    G4double WallTyvekInnerRadius;
    G4double WallTyvekOuterRadius;
    G4double WallTyvekBarrelLength;

    G4VisAttributes* MainWaterTankVis;
    G4Material* MainWaterTankMaterial;
    G4double MainWaterTankRadius;
    G4double MainWaterTankLength;

    G4VisAttributes* RockShellVis;
    G4Material* RockShellMaterial;
    G4double RockShellRadius;
    G4double RockShellLength;

    G4double NFrameCellsPerRow;
    G4double NRowsPerMainBlock;
    G4double NRowsPerBottomBlock;
    G4double NBlocksPerRing;


    void Print(){
      std::cout << "BlackTyvekInnerRadius = " << BlackTyvekInnerRadius << std::endl;
      std::cout << "BlackTyvekOuterRadius = " << BlackTyvekOuterRadius << std::endl;
      std::cout << "DeadSpaceInnerRadius = " << DeadSpaceInnerRadius << std::endl;
      std::cout << "DeadSpaceOuterRadius = " << DeadSpaceOuterRadius << std::endl;
      std::cout << "WhiteTyvekInnerRadius = " << WhiteTyvekInnerRadius << std::endl;
      std::cout << "WallTyvekOuterRadius = " << WallTyvekOuterRadius << std::endl;
      std::cout << "WallTyvekInnerRadius = " << WallTyvekInnerRadius << std::endl;
      std::cout << "WallTyvekOuterRadius = " << WallTyvekOuterRadius << std::endl;
      std::cout << "MainWaterTankRadius = " << MainWaterTankRadius << std::endl;
      std::cout << "MainWaterTankLength = " << MainWaterTankLength << std::endl;
      std::cout << "RockShellRadius = " << RockShellLength << std::endl;
      std::cout << "RockShellLength = " << RockShellLength << std::endl;
    }

};

RelasticPlacementConfiguration config;

G4LogicalVolume* BuildAndPlace_SinglePolyhedraTank(
  std::string name,
  G4double start_radius,
  G4double end_radius,
  G4double full_length,
  G4Material* material,
  G4ThreeVector& position,
  G4LogicalVolume* mother,
  G4VisAttributes* vis
)
{

  // G4Tubs* solid = new G4Tubs(name,
  //                             start_radius,
  //                             end_radius,
  //                             0.5*full_length,
  //                             0.*deg,
  //                             360.*deg);

  G4double zplane[2] = {-0.5*full_length,0.5*full_length};
  G4double rstart[2] = {start_radius,start_radius};
  G4double rend[2] = {end_radius,end_radius};

  G4Polyhedra* solid = new G4Polyhedra(name,
                                        0, // phi start
                                        360*deg,
                                        48*6, //NPhi-gon
                                        2,
                                        zplane,
                                        rstart,
                                        rend);

  G4LogicalVolume* logic = 
      new G4LogicalVolume(solid,
                          material,
                          name,
                          0,0,0);

  if (mother){
  G4PVPlacement* physical = new G4PVPlacement(0,
                          G4ThreeVector(0.,0.,0.),
                          logic,
                          name,
                          mother,
                          false,
                          0,
                          0); 
  }

  vis->SetForceLineSegmentsPerCircle(128);
  logic->SetVisAttributes(vis); 

  return logic;
}


G4LogicalVolume* WCSimDetectorConstruction::ConstructRealisticPlacement()
{
  G4cout << "**** Building Realistic HK Placement Detector ****" << G4endl;
   G4NistManager* nist = G4NistManager::Instance();

    // *************************
    // Legacy WCSim Code
    // *************************
    // WCSim has lots of calculations that fill DetectorConstruction Variables
    // To avoid having to search through the structure ALL variables needed 
    // are contained in the config structure.

    // New Configuration Calculation (eventually DB read)
    WCLength  = 70*m;

    WCIDRadius = 64.8*m/2;
    WCIDHeight = 65.751*m;
    WCODLateralWaterDepth    = 1.*m;
    WCODHeightWaterDepth     = 2.*m;
    WCODDeadSpace            = 600.*mm;
    WCODTyvekSheetThickness  = 1.*mm; // Quite standard I guess
    WCBlackSheetThickness = 2*mm;

    config.InnerDetectorVis = new G4VisAttributes(true, G4Colour(0.0,0.0,1.0,1.0));
    config.InnerDetectorVis->SetForceSolid(1);
    config.InnerDetectorMaterial = G4Material::GetMaterial("Water");
    config.InnerDetectorInnerRadius = 0;
    config.InnerDetectorOuterRadius = WCIDRadius;
    config.InnerDetectorBarrelLength = WCIDHeight;

    config.BlackTyvekVis = new G4VisAttributes(true, G4Colour(0.0,1.0,0.0,1.0));
    // config.BlackTyvekVis->SetForceWireframe(1);
    config.BlackTyvekVis->SetLineWidth(2);
    config.BlackTyvekVis->SetForceAuxEdgeVisible(0);
    config.BlackTyvekMaterial = G4Material::GetMaterial("Tyvek");
    config.BlackTyvekInnerRadius = WCIDRadius;
    config.BlackTyvekOuterRadius = WCIDRadius + WCBlackSheetThickness;
    config.BlackTyvekBarrelLength = WCIDHeight + 2*WCBlackSheetThickness;

    config.DeadSpaceVis = new G4VisAttributes(true, G4Colour(0.0,0.0,0.0,1.0));
    config.DeadSpaceVis->SetForceSolid(1);
    config.DeadSpaceMaterial = G4Material::GetMaterial("Water");
    config.DeadSpaceInnerRadius = config.BlackTyvekOuterRadius;
    config.DeadSpaceOuterRadius = config.BlackTyvekOuterRadius + WCODDeadSpace;
    config.DeadSpaceBarrelLength = config.BlackTyvekBarrelLength + 2*WCODDeadSpace;

    config.WhiteTyvekVis = new G4VisAttributes(true, G4Colour(1.0,1.0,1.0,1.0));
    config.WhiteTyvekVis->SetForceSolid(1);
    config.WhiteTyvekVis->SetLineWidth(2);
    config.WhiteTyvekVis->SetForceAuxEdgeVisible(0);
    config.WhiteTyvekMaterial = G4Material::GetMaterial("Tyvek");
    config.WhiteTyvekInnerRadius = config.DeadSpaceOuterRadius;
    config.WhiteTyvekOuterRadius = config.DeadSpaceOuterRadius + WCODTyvekSheetThickness;
    config.WhiteTyvekBarrelLength = config.DeadSpaceBarrelLength + 2*WCODTyvekSheetThickness;

    config.OuterDetectorVis = new G4VisAttributes(false, G4Colour(0.5,0.5,1.0,0.05));
    config.OuterDetectorMaterial = G4Material::GetMaterial("Water");
    config.OuterDetectorInnerRadius = config.WhiteTyvekOuterRadius;
    config.OuterDetectorOuterRadius = config.WhiteTyvekOuterRadius + WCODLateralWaterDepth;
    config.OuterDetectorBarrelLength = config.WhiteTyvekBarrelLength + 2*WCODHeightWaterDepth;

    config.WallTyvekVis = new G4VisAttributes(false, G4Colour(1.0,0.0,0.0,1.0));
    config.WallTyvekVis->SetForceWireframe(1);
    config.WallTyvekVis->SetLineWidth(0);
    config.WallTyvekVis->SetForceAuxEdgeVisible(0);
    config.WallTyvekMaterial = G4Material::GetMaterial("Tyvek");
    config.WallTyvekInnerRadius = config.OuterDetectorOuterRadius;
    config.WallTyvekOuterRadius = config.OuterDetectorOuterRadius + WCODTyvekSheetThickness;
    config.WallTyvekBarrelLength = config.OuterDetectorBarrelLength + 2*WCODTyvekSheetThickness;

    config.MainWaterTankVis = new G4VisAttributes(false, G4Colour(0.0,0.1,0.2,0.5));
    config.MainWaterTankRadius = config.WallTyvekOuterRadius+5*mm;
    config.MainWaterTankLength = config.WallTyvekBarrelLength+5*mm;
    config.MainWaterTankMaterial = G4Material::GetMaterial("Water");

    config.RockShellVis = new G4VisAttributes(false, G4Colour(0.5,0.5,1.0,0.1));
    config.RockShellMaterial = G4Material::GetMaterial("Water");
    config.RockShellRadius = config.MainWaterTankRadius + 30*cm;
    config.RockShellLength = config.MainWaterTankLength + 30*cm;

    config.NFrameCellsPerRow = 48;
    config.NRowsPerMainBlock = 16;
    config.NRowsPerBottomBlock = 12;
    config.NBlocksPerRing = 6;

    config.Print();

    // Our Hiierarchy is as follows
    G4ThreeVector CENTRAL_POS = G4ThreeVector();

    // 1. Rock (WC)
    // - OD WCBarrel 
    // - - White Tyvek Wall
    // - - - OD Water Space
    // - - - - OD PMTs
    // - - - - White Tyvek Frame
    // - - - - Black Tyvek Wall
    // - - - - - Dead Space
    // - - - - - - White Tyvek
    // - - - - - - - Water Inner
    // - - - - - - - - ID PMTs

    // 1. Rock (WC)
    // This is actually a rock solid, moved it instead of air.
    G4LogicalVolume* rockShellLogic = BuildAndPlace_SinglePolyhedraTank(
      "WC",
      0.0,
      config.RockShellRadius,
      config.RockShellLength,
      config.RockShellMaterial,
      CENTRAL_POS,
      NULL,
      config.RockShellVis
    );
 
    // 2. OD (WCBarrel)
    G4LogicalVolume* MainWaterTankLogic = BuildAndPlace_SinglePolyhedraTank(
      "WCBarrel",
      0.0,
      config.MainWaterTankRadius,
      config.MainWaterTankLength,
      config.MainWaterTankMaterial,
      CENTRAL_POS,
      rockShellLogic,
      config.MainWaterTankVis
    );

    // 3. White Tyvek Shell
    G4LogicalVolume* WallTyvekLogic = BuildAndPlace_SinglePolyhedraTank(
      "WallTyvek",
      0.0,
      config.WallTyvekOuterRadius,
      config.WallTyvekBarrelLength,
      config.WallTyvekMaterial,
      CENTRAL_POS,
      MainWaterTankLogic,
      config.WallTyvekVis
    );

    // 4. OuterDetectr
    G4LogicalVolume* OuterDetectorLogic = BuildAndPlace_SinglePolyhedraTank(
      "OuterDetector",
      0.0,
      config.OuterDetectorOuterRadius,
      config.OuterDetectorBarrelLength,
      config.OuterDetectorMaterial,
      CENTRAL_POS,
      WallTyvekLogic,
      config.OuterDetectorVis
    );

    // 5. WhiteTyvek
    G4LogicalVolume* WhiteTyvekLogic = BuildAndPlace_SinglePolyhedraTank(
      "WhiteTyvek",
      0.0,
      config.WhiteTyvekOuterRadius,
      config.WhiteTyvekBarrelLength,
      config.WhiteTyvekMaterial,
      CENTRAL_POS,
      OuterDetectorLogic,
      config.WhiteTyvekVis
    );

    // 6. Dead Space
    G4LogicalVolume* DeadSpaceLogic = BuildAndPlace_SinglePolyhedraTank(
      "DeadSpace",
      0.0,
      config.DeadSpaceOuterRadius,
      config.DeadSpaceBarrelLength,
      config.DeadSpaceMaterial,
      CENTRAL_POS,
      WhiteTyvekLogic,
      config.DeadSpaceVis
    );

    // 7. BlackTyvek
    G4LogicalVolume* BlackTyvekLogic = BuildAndPlace_SinglePolyhedraTank(
      "BlackTyvek",
      0.0,
      config.BlackTyvekOuterRadius,
      config.BlackTyvekBarrelLength,
      config.BlackTyvekMaterial,
      CENTRAL_POS,
      DeadSpaceLogic,
      config.BlackTyvekVis
    );

    // 8. InnerDetector
    G4LogicalVolume* InnerDetectorLogic = BuildAndPlace_SinglePolyhedraTank(
      "InnerDetector",
      0.0,
      config.InnerDetectorOuterRadius,
      config.InnerDetectorBarrelLength,
      config.InnerDetectorMaterial,
      CENTRAL_POS,
      BlackTyvekLogic,
      config.InnerDetectorVis
    );

    // Optionall inner phantom
    // G4LogicalVolume* InnerDetectorPhantomLogic = BuildAndPlace_SinglePolyhedraTank(
    //   "InnerDetectorPhantom",
    //   0.0,
    //   config.InnerDetectorOuterRadius-3*m,
    //   config.InnerDetectorBarrelLength-3*m,
    //   config.InnerDetectorMaterial,
    //   CENTRAL_POS,
    //   InnerDetectorLogic,
    //   config.InnerDetectorVis
    // );


    G4double twopi = 3.14159*2;
    G4double pi = 3.14159*2;

    // Temporary PMT Logics to place
    G4Material* pmt_mat = nist->FindOrBuildMaterial("G4_AIR");
    auto pmt20_solid = new G4Sphere("pmt20",0, 508*mm/2, 0, pi,0,pi);
    auto pmt20_logic = new G4LogicalVolume(pmt20_solid, pmt_mat, "pmt20");

    auto pmtMulti_solid = new G4Sphere("pmtMulti",0, 508*mm/2, 0, pi,0,pi);
    auto pmtMulti_logic = new G4LogicalVolume(pmtMulti_solid, pmt_mat, "pmtMulti");

    auto pmtod_solid = new G4Sphere("pmtod",0, 200*mm/2, 0, pi,0,pi);
    auto pmtod_logic = new G4LogicalVolume(pmtod_solid, pmt_mat, "pmtod");


    G4VisAttributes* showColor = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
    showColor->SetForceLineSegmentsPerCircle(12);
    showColor->SetForceAuxEdgeVisible(1);
    showColor->SetForceWireframe(1);
    showColor->SetLineWidth(3);


    G4VisAttributes* showColor2 = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
    showColor2->SetForceLineSegmentsPerCircle(12);
    showColor2->SetForceAuxEdgeVisible(1);
    showColor2->SetForceWireframe(1);
    showColor2->SetLineWidth(3);

    G4VisAttributes* showColor3 = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    showColor3->SetForceLineSegmentsPerCircle(12);
    showColor3->SetForceAuxEdgeVisible(1);
    showColor3->SetForceWireframe(1);
    showColor3->SetLineWidth(3);

    pmt20_logic->SetVisAttributes(showColor);
    pmtMulti_logic->SetVisAttributes(showColor2);
    pmtod_logic->SetVisAttributes(showColor3);


    // First make a single cell assembly, this will put the PMT flush against the wall
    G4double InnerTyvekRadius = config.InnerDetectorOuterRadius;

    auto pmt_cell_assembly = new G4AssemblyVolume();
    G4RotationMatrix* pmtcentralrotation = new G4RotationMatrix;
    pmtcentralrotation->rotateZ(180*deg);
    pmtcentralrotation->rotateY(270*deg);
    pmtcentralrotation->rotateX(90*deg);

    G4ThreeVector pmt_central_position = G4ThreeVector(InnerTyvekRadius,0.0,0.0);
    G4ThreeVector pmt_central_offset = G4ThreeVector(0.0,0.0,0.0); // -> Can be used for relative offset!
    pmt_cell_assembly->AddPlacedVolume(pmt20_logic, pmt_central_position, pmtcentralrotation);

    auto pmtmulti_cell_assembly = new G4AssemblyVolume();
    G4RotationMatrix* multipmtcentralrotation = new G4RotationMatrix;
    multipmtcentralrotation->rotateY(270*deg);
    multipmtcentralrotation->rotateX(270*deg);
    pmtmulti_cell_assembly->AddPlacedVolume(pmtMulti_logic, pmt_central_position, multipmtcentralrotation);

    auto odpmt_cell_assembly = new G4AssemblyVolume();
    G4ThreeVector odpmt_central_position = G4ThreeVector(config.WhiteTyvekOuterRadius+5*cm,0.0,0.0);
    G4RotationMatrix* odpmtcentralrotation = new G4RotationMatrix;
    odpmtcentralrotation->rotateY(270*deg);
    odpmt_cell_assembly->AddPlacedVolume(pmtod_logic, odpmt_central_position, odpmtcentralrotation);

    G4double RowSeperation = 1.4*m;

    int NSpacesInBlock = 48;
    int NSegments = 6*NSpacesInBlock;
    G4double phi_offset = twopi / NSegments;

    G4RotationMatrix* pmtsteprotation = new G4RotationMatrix;
    G4RotationMatrix* pmtsteprotation2 = new G4RotationMatrix;
    pmtsteprotation2->rotateZ(-phi_offset/2);
    pmtsteprotation->rotateZ(phi_offset/2);//align the PMT with the Cell

    G4ThreeVector poslow = G4ThreeVector(0,0.0,-RowSeperation/4);
    G4ThreeVector posup = G4ThreeVector(0,0.0,RowSeperation/4);

    auto frame_block_assembly = new G4AssemblyVolume();
    frame_block_assembly->AddPlacedAssembly(pmt_cell_assembly, posup, pmtsteprotation  );
    frame_block_assembly->AddPlacedAssembly(pmt_cell_assembly, poslow, pmtsteprotation2  );

    auto frame_block_assembly_withpmt = new G4AssemblyVolume();  
    frame_block_assembly_withpmt->AddPlacedAssembly(pmt_cell_assembly, posup, pmtsteprotation  );
    frame_block_assembly_withpmt->AddPlacedAssembly(pmt_cell_assembly, poslow, pmtsteprotation2  );
    frame_block_assembly_withpmt->AddPlacedAssembly(pmtmulti_cell_assembly, poslow, pmtsteprotation  );

    auto frame_block_assembly_odpmt = new G4AssemblyVolume();
    G4RotationMatrix* odsteprotation = new G4RotationMatrix;
    G4ThreeVector posod = G4ThreeVector(0,0.0,-RowSeperation/2);
    frame_block_assembly_odpmt->AddPlacedAssembly(odpmt_cell_assembly, posod, odsteprotation);

    // Make the Three types of block rows
    auto block_row_nomultipmts = new G4AssemblyVolume();
    auto block_row_index1multipmts = new G4AssemblyVolume();
    auto block_row_index4multipmts = new G4AssemblyVolume();

    auto block_row_odzeroindex = new G4AssemblyVolume();
    auto block_row_odoffsetindex = new G4AssemblyVolume();
    auto block_row_odfourindex = new G4AssemblyVolume();


    G4ThreeVector poscentral;
    auto pmtcentralrotation2 = new G4RotationMatrix();
    for (int i = 0; i < NSpacesInBlock/2; i++){
      pmtcentralrotation2->rotateZ(phi_offset*2);//align the PMT with the Cell
      
      block_row_nomultipmts->AddPlacedAssembly(frame_block_assembly, poscentral, pmtcentralrotation2);

      if ((i-1) % 6 != 0){
        block_row_index1multipmts->AddPlacedAssembly(frame_block_assembly, poscentral, pmtcentralrotation2);
      } else {
        block_row_index1multipmts->AddPlacedAssembly(frame_block_assembly_withpmt, poscentral, pmtcentralrotation2);
      }

      if ((i-4) % 6 != 0){
        block_row_index4multipmts->AddPlacedAssembly(frame_block_assembly, poscentral, pmtcentralrotation2);
      } else {
        block_row_index4multipmts->AddPlacedAssembly(frame_block_assembly_withpmt, poscentral, pmtcentralrotation2);
      }

      if (i % 2 == 0) block_row_odzeroindex->AddPlacedAssembly(frame_block_assembly_odpmt, poscentral, pmtcentralrotation2);
      if ((i-3) % 4 == 0) block_row_odoffsetindex->AddPlacedAssembly(frame_block_assembly_odpmt, poscentral, pmtcentralrotation2);
      if ((i-1) % 4 == 0) block_row_odfourindex->AddPlacedAssembly(frame_block_assembly_odpmt, poscentral, pmtcentralrotation2);

   }

    // Make the barrel block 
    auto block_assembly = new G4AssemblyVolume();
    auto block_bottom = new G4AssemblyVolume();

    auto odblock_assembly = new G4AssemblyVolume();
    auto odblock_bottom = new G4AssemblyVolume();

    G4RotationMatrix* blockrotation = new G4RotationMatrix();


    // Hard coded for now
    G4ThreeVector block_pos = G4ThreeVector(0.0,0.0,-0.5*RowSeperation);
    block_assembly->AddPlacedAssembly(block_row_nomultipmts, block_pos, blockrotation);

    block_pos = G4ThreeVector(0.0,0.0,-1.5*RowSeperation);
    block_assembly->AddPlacedAssembly(block_row_index1multipmts, block_pos, blockrotation);

    block_pos = G4ThreeVector(0.0,0.0,-2.5*RowSeperation);
    block_assembly->AddPlacedAssembly(block_row_nomultipmts, block_pos, blockrotation);

    block_pos = G4ThreeVector(0.0,0.0,-3.5*RowSeperation);
    block_assembly->AddPlacedAssembly(block_row_index4multipmts, block_pos, blockrotation);

    block_pos = G4ThreeVector(0.0,0.0,-4.5*RowSeperation);
    block_assembly->AddPlacedAssembly(block_row_nomultipmts, block_pos, blockrotation);

    block_pos = G4ThreeVector(0.0,0.0,-5.5*RowSeperation);
    block_assembly->AddPlacedAssembly(block_row_index1multipmts, block_pos, blockrotation);

    block_pos = G4ThreeVector(0.0,0.0,-6.5*RowSeperation);
    block_assembly->AddPlacedAssembly(block_row_nomultipmts, block_pos, blockrotation);

    block_pos = G4ThreeVector(0.0,0.0,-7.5*RowSeperation);
    block_assembly->AddPlacedAssembly(block_row_index4multipmts, block_pos, blockrotation);


    // Hard coded for now
    block_pos = G4ThreeVector(0.0,0.0,-0.5*RowSeperation);
    block_bottom->AddPlacedAssembly(block_row_nomultipmts, block_pos, blockrotation);

    block_pos = G4ThreeVector(0.0,0.0,-1.5*RowSeperation);
    block_bottom->AddPlacedAssembly(block_row_index1multipmts, block_pos, blockrotation);

    block_pos = G4ThreeVector(0.0,0.0,-2.5*RowSeperation);
    block_bottom->AddPlacedAssembly(block_row_nomultipmts, block_pos, blockrotation);

    block_pos = G4ThreeVector(0.0,0.0,-3.5*RowSeperation);
    block_bottom->AddPlacedAssembly(block_row_index4multipmts, block_pos, blockrotation);

    block_pos = G4ThreeVector(0.0,0.0,-4.5*RowSeperation);
    block_bottom->AddPlacedAssembly(block_row_nomultipmts, block_pos, blockrotation);

    block_pos = G4ThreeVector(0.0,0.0,-5.5*RowSeperation);
    block_bottom->AddPlacedAssembly(block_row_index1multipmts, block_pos, blockrotation);

    // OD Block
    G4ThreeVector odblock_pos = G4ThreeVector(0.0,0.0,-0.5*RowSeperation);
    odblock_assembly->AddPlacedAssembly(block_row_odzeroindex, odblock_pos, blockrotation);

    odblock_pos = G4ThreeVector(0.0,0.0,-1.5*RowSeperation);
    odblock_assembly->AddPlacedAssembly(block_row_odoffsetindex, odblock_pos, blockrotation);

    odblock_pos = G4ThreeVector(0.0,0.0,-2.5*RowSeperation);
    odblock_assembly->AddPlacedAssembly(block_row_odzeroindex, odblock_pos, blockrotation);

    odblock_pos = G4ThreeVector(0.0,0.0,-3.5*RowSeperation);
    odblock_assembly->AddPlacedAssembly(block_row_odfourindex, odblock_pos, blockrotation);

    odblock_pos = G4ThreeVector(0.0,0.0,-4.5*RowSeperation);
    odblock_assembly->AddPlacedAssembly(block_row_odzeroindex, odblock_pos, blockrotation);

    odblock_pos = G4ThreeVector(0.0,0.0,-5.5*RowSeperation);
    odblock_assembly->AddPlacedAssembly(block_row_odoffsetindex, odblock_pos, blockrotation);

    odblock_pos = G4ThreeVector(0.0,0.0,-6.5*RowSeperation);
    odblock_assembly->AddPlacedAssembly(block_row_odzeroindex, odblock_pos, blockrotation);

    odblock_pos = G4ThreeVector(0.0,0.0,-7.5*RowSeperation);
    odblock_assembly->AddPlacedAssembly(block_row_odfourindex, odblock_pos, blockrotation);

     // OD bottom Block
    odblock_pos = G4ThreeVector(0.0,0.0,-0.5*RowSeperation);
    odblock_bottom->AddPlacedAssembly(block_row_odzeroindex, odblock_pos, blockrotation);

    odblock_pos = G4ThreeVector(0.0,0.0,-1.5*RowSeperation);
    odblock_bottom->AddPlacedAssembly(block_row_odoffsetindex, odblock_pos, blockrotation);

    odblock_pos = G4ThreeVector(0.0,0.0,-2.5*RowSeperation);
    odblock_bottom->AddPlacedAssembly(block_row_odzeroindex, odblock_pos, blockrotation);

    odblock_pos = G4ThreeVector(0.0,0.0,-3.5*RowSeperation);
    odblock_bottom->AddPlacedAssembly(block_row_odfourindex, odblock_pos, blockrotation);

    odblock_pos = G4ThreeVector(0.0,0.0,-4.5*RowSeperation);
    odblock_bottom->AddPlacedAssembly(block_row_odoffsetindex, odblock_pos, blockrotation);

    odblock_pos = G4ThreeVector(0.0,0.0,-5.5*RowSeperation);
    odblock_bottom->AddPlacedAssembly(block_row_odfourindex, odblock_pos, blockrotation);



    G4double imprint_spacing = 8*RowSeperation;

    // Inprint Columns first then rings
    for (int i = 0; i < 5; i++){
      G4RotationMatrix* imprint_rot = new G4RotationMatrix();
      G4ThreeVector imprint_pos = G4ThreeVector(0.0,0.0,-i*imprint_spacing + WCIDHeight/2-60*cm);
      for (int j = 0; j < 6; j++){
        imprint_rot->rotateZ( twopi / 6 );
        block_assembly->MakeImprint(InnerDetectorLogic, imprint_pos, imprint_rot);
        odblock_assembly->MakeImprint(OuterDetectorLogic, imprint_pos, imprint_rot);
      }
    }

    G4RotationMatrix* imprint_rot2 = new G4RotationMatrix();
    G4ThreeVector imprint_pos2 = G4ThreeVector(0.0,0.0,-5*imprint_spacing + WCIDHeight/2-60*cm);
    for (int j = 0; j < 6; j++){
      imprint_rot2->rotateZ( twopi / 6 );
      block_bottom->MakeImprint(InnerDetectorLogic, imprint_pos2, imprint_rot2);
      odblock_bottom->MakeImprint(OuterDetectorLogic, imprint_pos2, imprint_rot2);
    }

    // End Cap Assembly
    G4AssemblyVolume* endcap_assembly = new G4AssemblyVolume();
    G4AssemblyVolume* endcap_assembly_od = new G4AssemblyVolume();

    G4double cap_offset_x = 0.52*m;
    G4double cap_offset_y = 0.52*m;
    G4RotationMatrix* endcapmtprotation = new G4RotationMatrix();
    endcapmtprotation->rotateX(180*deg);

    int nplaced = 0;
    for (int i = -62; i < 63; i++ ){
      for (int j = -62; j < 63; j++ ){
        G4ThreeVector pmtpos = G4ThreeVector( i*cap_offset_x, j*cap_offset_y, 0.0 );
        if (i % 2 == 0) continue;
        if (j % 2 == 0) continue;
        if (pmtpos.mag() > WCIDRadius - 40*cm) continue;
        nplaced++;
        endcap_assembly->AddPlacedVolume(pmt20_logic, pmtpos, endcapmtprotation);
      }
    }

    // OD Cap
    for (int i = -62; i < 63; i++ ){
      for (int j = -62; j < 63; j++ ){
        G4ThreeVector pmtpos = G4ThreeVector( i*cap_offset_x, j*cap_offset_y, 0.0 );
        if (!(i % 2 == 0)) continue;
        if (!(j % 3 == 0)) continue;
        // if (i % 4 == 0 && (j % 6 == 0)) continue;
        bool valid = false;
        if ((i % 2 == 0 && ((j-3) % 6 == 0) && ((i) % 4 == 0))) valid = true;
        if ((j % 3 == 0 && ((j) % 6 == 0) && ((i-2) % 4 == 0))) valid = true;
        if (!valid) continue;
        if (pmtpos.mag() > WCIDRadius - 40*cm) continue;
        endcap_assembly_od->AddPlacedVolume(pmtod_logic, pmtpos, endcapmtprotation);
      }
    }


    nplaced = 0;
    for (int i = -62; i < 63; i++ ){
      for (int j = -62; j < 63; j++ ){
        G4ThreeVector pmtpos = G4ThreeVector( i*cap_offset_x, j*cap_offset_y, 0.0 );
        if (i % 2 == 1) continue;
        if (j % 2 == 1) continue;
        if (!(i % 6 == 0 && j % 6 == 0 && !(i % 12 == 0 && j % 12 == 0))) continue;
        if (pmtpos.mag() > WCIDRadius - 40*cm) continue;
        nplaced++;
        endcap_assembly->AddPlacedVolume(pmtMulti_logic, pmtpos, endcapmtprotation);
      }
    }

    G4ThreeVector topcappos = G4ThreeVector(0.0,0.0,WCIDHeight/2);
    G4RotationMatrix* topcaprot = new G4RotationMatrix();
    endcap_assembly->MakeImprint( InnerDetectorLogic, topcappos, topcaprot);

    G4ThreeVector botcappos = G4ThreeVector(0.0,0.0,-WCIDHeight/2);
    G4RotationMatrix* botcaprot = new G4RotationMatrix();
    botcaprot->rotateX(180*deg);
    endcap_assembly->MakeImprint( InnerDetectorLogic, botcappos, botcaprot);

    G4ThreeVector topcapposod = G4ThreeVector(0.0,0.0,config.WhiteTyvekBarrelLength/2);
    G4RotationMatrix* topcaprotod = new G4RotationMatrix();
    endcap_assembly_od->MakeImprint( OuterDetectorLogic, topcapposod, topcaprotod);

    G4ThreeVector botcapposod = G4ThreeVector(0.0,0.0,-config.WhiteTyvekBarrelLength/2);
    G4RotationMatrix* botcaprotod = new G4RotationMatrix();
    botcaprotod->rotateX(180*deg);
    endcap_assembly_od->MakeImprint( OuterDetectorLogic, botcapposod, botcaprotod);


    // -------------------------------------
    // ID PMT Creation (copoied from WCSIM)
    // -------------------------------------
    std::cout << "Inner detector has : " << InnerDetectorLogic->GetNoDaughters() << std::endl;
    std::cout << "Outer detector has : " << OuterDetectorLogic->GetNoDaughters() << std::endl;

    G4LogicalVolume* logicWCPMT;
    if(nID_PMTs<=1) logicWCPMT = ConstructPMT(WCPMTName, WCIDCollectionName,"tank",nID_PMTs);
    else logicWCPMT = ConstructMultiPMT(WCPMTName, WCIDCollectionName,"tank",nID_PMTs);
    if(!logicWCPMT){
      G4cerr << "Overlapping PMTs in multiPMT" << G4endl;
      return NULL; 
    }
    G4LogicalVolume* logicWCPMT2 = nullptr;
  #ifdef DEBUG
    G4cout << "HYBRID2 = " << hybrid << G4endl;
  #endif
    if(hybrid){
      G4cout<<"First type of PMT LV is created. Now creating the LV for the second type of PMT"<<G4endl;
      if(nID_PMTs2<=1) logicWCPMT2 = ConstructPMT(WCPMTName2, WCIDCollectionName2,"tankPMT2",nID_PMTs2);
      else logicWCPMT2 = ConstructMultiPMT(WCPMTName2, WCIDCollectionName2,"tankPMT2",nID_PMTs2);
      if(!logicWCPMT2){
        G4cerr << "Overlapping PMTs in multiPMT" << G4endl;
        return NULL; 
      }
    }

    G4String pmtname = "WCMultiPMT";


    // -------------------------------------
    // ID PMT Placement
    // -------------------------------------
    int ndaughters = InnerDetectorLogic->GetNoDaughters();
    int copyno = 1;
    std::vector<G4Transform3D> positions;

    int removed = 0;
    for (int i = 0; i < ndaughters; i++){
      G4VPhysicalVolume* vol = InnerDetectorLogic->GetDaughter(i-removed);
      G4Transform3D aTransform = G4Transform3D(*(vol->GetObjectRotation()), vol->GetObjectTranslation());

      if ( vol->GetLogicalVolume() == pmt20_logic ) {

        new G4PVPlacement(
            aTransform,
            logicWCPMT,              // pmt20
            pmtname,       // its name
            InnerDetectorLogic,       // its mother volume
            false,                    // no boolean operations
            copyno++,
            false
        );   
         
        InnerDetectorLogic->RemoveDaughter(vol);
        delete vol;
        removed++;

      } else if ( vol->GetLogicalVolume() == pmtMulti_logic ) {

          new G4PVPlacement(
              aTransform,
              logicWCPMT2,              // pmt20
              pmtname,       // its name
              InnerDetectorLogic,       // its mother volume
              false,                    // no boolean operations
              copyno++,
              false
          );   

          logicWCPMT2->SetVisAttributes(showColor2);          
          InnerDetectorLogic->RemoveDaughter(vol);
          delete vol;
          removed++;
      }
    }

    // -------------------------------------
    // OD PMT Creation and Placement
    // -------------------------------------
    logicWCODWLSAndPMT = ConstructPMTAndWLSPlate(WCPMTODName, WCODCollectionName, "OD");

    removed = 0;
    ndaughters = OuterDetectorLogic->GetNoDaughters();
    for (int i = 0; i < ndaughters; i++){
      G4VPhysicalVolume* vol = OuterDetectorLogic->GetDaughter(i-removed);
      G4Transform3D aTransform = G4Transform3D(*(vol->GetObjectRotation()), vol->GetObjectTranslation());

      if ( vol->GetLogicalVolume() == pmtod_logic ) {

        // if (copyno < 1000){
        new G4PVPlacement(
            aTransform,
            logicWCODWLSAndPMT,              // pmt20
            "WCBorderCellODContainer",       // its name
            OuterDetectorLogic,       // its mother volume
            false,                    // no boolean operations
            copyno++,
            false
        );   
                
        OuterDetectorLogic->RemoveDaughter(vol);
        delete vol;
        removed++;

      }

    }
  
    // Optical Surfaces
    // Main Optical Skin Surfaces are handled inside PMTs themselves.
    
    // Nested structure means we should just need one skinsurfaces for all the tyveks.
    new G4LogicalSkinSurface("WallTyvekSurface",WallTyvekLogic,OpWaterTySurface);
    new G4LogicalSkinSurface("WhiteTyvekSurface",WhiteTyvekLogic,OpWaterTySurface);
    new G4LogicalSkinSurface("BlackTyvekSurface",BlackTyvekLogic,OpWaterBSSurface);

    std::cout << "FINISHED Realistic Placement" << std::endl;
    std::cout << "FINISHED Realistic Placement" << std::endl;
    std::cout << "FINISHED Realistic Placement" << std::endl;
    std::cout << "FINISHED Realistic Placement" << std::endl;
    std::cout << "FINISHED Realistic Placement" << std::endl;
    std::cout << "FINISHED Realistic Placement" << std::endl;
    std::cout << "FINISHED Realistic Placement" << std::endl;
    std::cout << "FINISHED Realistic Placement" << std::endl;
    std::cout << "FINISHED Realistic Placement" << std::endl;

    return rockShellLogic;
}

