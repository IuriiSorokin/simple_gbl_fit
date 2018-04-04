
#include "RKTrackRep.h"
#include "GblFitter.h"
#include "TGeoElement.h"
#include "TGeoManager.h"
#include "TGeoMaterialInterface.h"
#include "MaterialEffects.h"
#include "FieldManager.h"
#include "ConstField.h"
#include "PlanarMeasurement.h"
#include "MeasuredStateOnPlane.h"
#include "TList.h"
#include "TRandom.h"

using namespace genfit;


int main()
{
    // ===========================================================================
    //     Settings
    // ===========================================================================

    constexpr size_t N = 5; // number of detector planes

    std::array<SharedPlanePtr, N> planes = { std::make_shared<DetPlane>( TVector3(0,0,10), TVector3( 1, 0, 0), TVector3(0, 1, 0) ),
                                             std::make_shared<DetPlane>( TVector3(0,0,30), TVector3(-1, 0, 0), TVector3(0, 1, 0) ),   // <-- notice the plane orientation
                                             std::make_shared<DetPlane>( TVector3(0,0,50), TVector3( 1, 0, 0), TVector3(0, 1, 0) ),
                                             std::make_shared<DetPlane>( TVector3(0,0,70), TVector3( 1, 0, 0), TVector3(0, 1, 0) ),
                                             std::make_shared<DetPlane>( TVector3(0,0,90), TVector3( 1, 0, 0), TVector3(0, 1, 0) ) };

    TVector3 true_position = { 0, 0, 0}; // true initial position of the
    TVector3 seed_position = { 0, 0, 0}; // seed for the fitting

    TVector3 true_momentum = { -0.01, 0,  0.1  };
    TVector3 seed_momentum = {     0, 0,  0.2 };

    TVector3 field = {0, 1, 0};

    // ===========================================================================
    //     Construct the world
    // ===========================================================================

    std::cout << " === Constructing the world ===" << std::endl;

    const double world_halfsize = 100;
    const double plane_thickness = 0.1;

    new TGeoManager("GenfitGeometry", "GenfitGeometry");

    TGeoManager::SetVerboseLevel(-1);

    genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

    TGeoElement * siliconElement = gGeoManager->GetElementTable()->FindElement("Silicon");
    TGeoElement * vacuumElement = gGeoManager->GetElementTable()->FindElement("Vacuum");

    Double_t mediumParameters[10];
    mediumParameters[0]=0.;//sensitive volume flag
    mediumParameters[1]=1.;//magnetic field flag
    mediumParameters[2]=30.;//max fiel in kGauss
    mediumParameters[3]=0.1;//maximal angular dev. due to field
    mediumParameters[4]=0.01;//max step allowed (in cm)
    mediumParameters[5]=1.e-5;//max fractional energy loss
    mediumParameters[6]=1.e-3;//boundary crossing accuracy
    mediumParameters[7]=1.e-5;//minimum step
    mediumParameters[8]=0.;//not defined
    mediumParameters[9]=0.;//not defined

    Double_t siliconDensity_g_cm3 = 2.3290;
    TGeoMaterial * siliconMaterial = new TGeoMaterial("siliconMaterial", siliconElement, siliconDensity_g_cm3);
    siliconMaterial->SetRadLen( 7777 ); // this should force ROOT to recalculate the radiation and interaction length. 7777 is an arbitrary positive number.
    std::cout << "Silicon radiation length = " << siliconMaterial->GetRadLen() << std::endl;
    TGeoMedium * siliconMedium = new TGeoMedium("siliconMedium", gGeoManager->GetListOfMedia()->GetSize(), siliconMaterial, mediumParameters);

    Double_t vacuumDensity_g_cm3 = 1e-27;
    TGeoMaterial * vacuumMaterial = new TGeoMaterial("vacuumMaterial", vacuumElement, vacuumDensity_g_cm3);
    vacuumMaterial->SetRadLen( 7777 ); // this should force ROOT to recalculate the radiation and interaction length. 7777 is an arbitrary positive number.
    std::cout << "Vacuum radiation length = " << vacuumMaterial->GetRadLen() << std::endl;
    TGeoMedium * vacuumMedium = new TGeoMedium("vacuumMedium", gGeoManager->GetListOfMedia()->GetSize(), vacuumMaterial, mediumParameters);

    // ===== Physical volumes =====
    TGeoVolume * topVolume = gGeoManager->MakeBox("TopVolume", vacuumMedium, world_halfsize, world_halfsize, world_halfsize);
    gGeoManager->SetTopVolume(topVolume);

    for( size_t i_plane = 0; i_plane < N; ++ i_plane ) {
        TGeoVolume * siliconPlane = gGeoManager->MakeBox("SiliconPlane", siliconMedium, world_halfsize, world_halfsize, plane_thickness / 2 );
        const auto& pos = planes.at(i_plane)->getO();
        topVolume->AddNode(siliconPlane, 0, new TGeoTranslation( pos.x(), pos.y(), pos.z() ) );
    }

    gGeoManager->CloseGeometry();

    genfit::FieldManager::getInstance()->init( new genfit::ConstField( field ) );

    // ===========================================================================
    //     Construct the measurements
    // ===========================================================================

    std::cout << " === Constructing the measurements ===" << std::endl;

    std::array< std::unique_ptr<PlanarMeasurement>, N> measurements;

    const int  electron_pdg_id = 11;

    RKTrackRep rep( electron_pdg_id ); // electron, propagation direction -- auto
    MeasuredStateOnPlane state( &rep );
    state.setPosMom( true_position, true_momentum );

    for( size_t i_plane = 0; i_plane < N; ++i_plane ) {
        std::cout << "Extrapolating to plane " << i_plane << std::endl;
        std::cout << "Init pos: "; state.getPos().Print();
        std::cout << "Init mom: "; state.getMom().Print();
        try {
            state.extrapolateToPlane( planes.at(i_plane), false, true );
        }
        catch( Exception& e ) {
            std::cout << "Extrapolation failed. Genfit exception:" << std::endl;
            e.info();
        }

        std::cout << "Ext. pos: "; state.getPos().Print();
        std::cout << "Ext. mom: "; state.getMom().Print();

        TMatrixDSym pos_cov(2);
        pos_cov(0,0) = state.getCov()(3,3);
        pos_cov(0,1) = state.getCov()(3,4);
        pos_cov(1,0) = state.getCov()(4,3);
        pos_cov(1,1) = state.getCov()(4,4);

        TVectorD measured_pos(2);
        measured_pos(0) = state.getPos().X() + 0.1 * gRandom->Gaus( 0, sqrt(pos_cov(0,0)) );
        measured_pos(1) = state.getPos().Y() + 0.1 * gRandom->Gaus( 0, sqrt(pos_cov(1,1)) );

        std::cout << "Meas.pos: "; TVector3(measured_pos(0), measured_pos(1), planes.at(i_plane)->getO().z() ).Print();
        std::cout << "Meas.err: "; ( TVector3(measured_pos(0), measured_pos(1), planes.at(i_plane)->getO().z() ) - state.getPos() ).Print();

        measurements.at(i_plane).reset( new PlanarMeasurement(measured_pos,pos_cov, i_plane, i_plane, nullptr ) );
        measurements.at(i_plane)->setPlane( planes.at(i_plane) );
    }

    // ===========================================================================
    //     Fitting
    // ===========================================================================

    std::cout << " === Fitting === " << std::endl;

    Track track( new genfit::RKTrackRep( electron_pdg_id ), seed_position, seed_momentum );

    for( size_t i_plane = 0; i_plane < N; ++i_plane ) {
        track.insertMeasurement( measurements.at(i_plane).release(), i_plane );
    }

    GblFitter fitter;

    fitter.processTrack( &track, false );

    //    // More fitting iterations
    //    for( int i = 0; i < 2; ++i ) {
    //        track.setStateSeed( track.getFittedState(0).get6DState() );
    //        fitter.processTrack( &track, false );
    //    }

    for( unsigned int iState = 0; iState < track.getNumPoints(); ++iState ) {
        std::cout << "State " << iState << std::endl;
        std::cout << "Fit. pos: "; track.getFittedState(iState).getPos().Print();
        std::cout << "Fit. mom: "; track.getFittedState(iState).getMom().Print();
    }

    return 0;
}


