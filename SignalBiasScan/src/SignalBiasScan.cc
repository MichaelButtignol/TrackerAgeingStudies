// -*- C++ -*-
//
// Package:    SignalBiasScan
// Class:      SignalBiasScan
// 
/**\class SignalBiasScan SignalBiasScan.cc SignalBiasScan/src/SignalBiasScan.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jeremy Andrea
//         Created:  Fri Apr 16 16:35:57 CEST 2010
// $Id: SignalBiasScan.cc,v 1.7 2013/02/19 12:53:50 jlagram Exp $
//      SignalBiasScan.cc,v 2.0 2012/09/03 J-L Agram
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"



#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include <Geometry/CommonTopologies/interface/Topology.h>
#include <Geometry/CommonTopologies/interface/StripTopology.h>

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include <DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h>
#include <DataFormats/SiStripDetId/interface/SiStripDetId.h>
#include <DataFormats/SiStripDetId/interface/TIBDetId.h>
#include <DataFormats/SiStripDetId/interface/TIDDetId.h>
#include <DataFormats/SiStripDetId/interface/TOBDetId.h>
#include <DataFormats/SiStripDetId/interface/TECDetId.h>
#include "DataFormats/DetId/interface/DetId.h"
#include <DataFormats/SiPixelCluster/interface/SiPixelCluster.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>


#include <DataFormats/TrackerRecHit2D/interface/SiTrackerMultiRecHit.h>
#include <DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h>

#include <RecoTracker/TransientTrackingRecHit/interface/TSiStripRecHit2DLocalPos.h>
#include <RecoTracker/TransientTrackingRecHit/interface/TSiTrackerMultiRecHit.h>
#include <RecoTracker/TransientTrackingRecHit/interface/TSiPixelRecHit.h>
#include <RecoTracker/TransientTrackingRecHit/interface/TSiStripRecHit1D.h>
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include <TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h>
#include <TrackingTools/PatternTools/interface/Trajectory.h>
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include <RecoLocalTracker/SiStripClusterizer/interface/SiStripClusterInfo.h>
#include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h" 


#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"


#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"


#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TTree.h"

#include "../interface/TreeEvent.h"

//
// function
//

bool isUsedInSmallBiasScan ( int detid ) {
  
  // TIB - 1.4.2.5, until run ??
  if(detid==369121605 || 
     detid==369121606 || 
     detid==369121614 || 
     detid==369121613 || 
     detid==369121610 ||
     detid==369121609 ) return true;
	 
  // TIB - 1.2.2.1, from run 193541
  if(detid==369121390 || 
     detid==369121382 || 
     detid==369121386 || 
     detid==369121385 || 
     detid==369121389 ||
     detid==369121381 ) return true;
  
  // TIB + 1.6.2.5
  if(detid==369125870 || 
     detid==369125862 || 
     detid==369125866 || 
     detid==369125869 || 
     detid==369125865 || // module not powered and masked
     detid==369125861 ) return true;
  
  // TOB + 4.3.3.8
  if(detid==436281512 || 
     detid==436281528 || 
     detid==436281508 || 
     detid==436281524 || 
     detid==436281520 || 
     detid==436281516 ) return true;
  
  // TEC
  if(detid==470148196 || 
     detid==470148200 || 
     detid==470148204 || 
     detid==470148228 || 
     detid==470148232 || 
     detid==470148236 || 
     detid==470148240 || 
     detid==470148261 || 
     detid==470148262 || 
     detid==470148265 || 
     detid==470148266 || 
     detid==470148292 || 
     detid==470148296 || 
     detid==470148300 || 
     detid==470148304 || 
     detid==470148324 || 
     detid==470148328 || 
     detid==470148332 || 
     detid==470148336 || 
     detid==470148340 ) return true; 


  return false;
}


//
// class declaration
//

class SignalBiasScan : public edm::EDAnalyzer {
   public:
      explicit SignalBiasScan(const edm::ParameterSet&);
      ~SignalBiasScan();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
  

      // ----------member data ---------------------------
  edm::InputTag trackLabel_;
  edm::ESHandle<TrackerGeometry> tracker_;
  edm::InputTag tkTraj_;
  edm::InputTag labelTrajToTrack_;
  edm::InputTag siStripClusters_;
  edm::InputTag primaryVertexColl_;
  
  bool fullHitInfo_TIB_;
  bool fullHitInfo_TOB_;
  bool fullHitInfo_TID_;
  bool fullHitInfo_TEC_;
  
  bool isSmallBiasScan_;


  edm::ESHandle<TransientTrackBuilder> theTTrackBuilder;
  //edm::ESHandle<GeometricSearchTracker> theTracker;
  edm::ESHandle<MagneticField> bField;  
  edm::ESHandle<Propagator> thePropagator;

  TTree *smalltree;
  
  edm::Service<TFileService> fs;

  TreeEvent *treeevent;
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SignalBiasScan::SignalBiasScan(const edm::ParameterSet& iConfig) {

  //now do what ever initialization is needed
  trackLabel_            = iConfig.getParameter<edm::InputTag>     ("trackLabel"      );
  tkTraj_                = iConfig.getParameter<edm::InputTag>     ("tkTraj"          );
  labelTrajToTrack_      = iConfig.getParameter<edm::InputTag>     ("labelTrajToTrack");
  siStripClusters_       = iConfig.getParameter<edm::InputTag>     ("siStripClusters" );
  primaryVertexColl_     = iConfig.getParameter<edm::InputTag>   ("primaryVertexColl"  );
  
  fullHitInfo_TIB_		= iConfig.getParameter< bool >   ( "fullHitInfo_TIB" );
  fullHitInfo_TOB_		= iConfig.getParameter< bool >   ( "fullHitInfo_TOB" );
  fullHitInfo_TID_		= iConfig.getParameter< bool >   ( "fullHitInfo_TID" );
  fullHitInfo_TEC_		= iConfig.getParameter< bool >   ( "fullHitInfo_TEC" );
  
  isSmallBiasScan_ 		= iConfig.getUntrackedParameter< bool > 	( "isSmallBiasScan", false );
  
  smalltree = fs->make<TTree>("ttree", "ttree");
  
  treeevent = new TreeEvent();
  
  smalltree->Branch("event", "TreeEvent", &treeevent, 32000, 3);
  
}


SignalBiasScan::~SignalBiasScan() {
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SignalBiasScan::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  
  //------------------
  //get primary vertex
  //------------------
  edm::Handle<reco::VertexCollection> primaryVertex;
  iEvent.getByLabel(primaryVertexColl_,primaryVertex);
  
  // Check there is at least one good primary vertex
  bool foundGoodVertex = true; 
  for(unsigned int i = 0; i< primaryVertex->size() ; i++){    
    if(!(*primaryVertex)[i].isFake() && (*primaryVertex)[i].tracksSize()>2 
       && fabs( (*primaryVertex)[i].z()) <= 24 && fabs( (*primaryVertex)[i].position().rho()) <= 2) foundGoodVertex = true; 
  }
  
  
  
  if(foundGoodVertex){

    edm::Handle<edm::View<reco::Track> > recoTrackHandle;
    iEvent.getByLabel(trackLabel_, recoTrackHandle);
    const edm::View<reco::Track> & recoTrackCollection = *(recoTrackHandle.product()); 
    edm::Handle<std::vector<Trajectory> > TrajectoryCollection;
    iEvent.getByLabel(tkTraj_,TrajectoryCollection);
    edm::Handle<TrajTrackAssociationCollection> trajTrackAssociationHandle;
    iEvent.getByLabel(labelTrajToTrack_,trajTrackAssociationHandle);
    
    edm::ESHandle<TrackerGeometry> theTrackerGeometry;
    iSetup.get<TrackerDigiGeometryRecord>().get( theTrackerGeometry );  
   
    edm::ESHandle<StripClusterParameterEstimator> parameterestimator;
    iSetup.get<TkStripCPERecord>().get("StripCPEfromTrackAngle", parameterestimator); 
    const StripClusterParameterEstimator &stripcpe(*parameterestimator);
    
	// Set event infos
    treeevent->reset();
    treeevent->run_nr = iEvent.id().run();
    treeevent->ev_nr  = iEvent.id().event();
    treeevent->ev_timestamp  = iEvent.time().unixTime();      

    // Loop on tracks
    for(TrajTrackAssociationCollection::const_iterator it = trajTrackAssociationHandle->begin(); it!=trajTrackAssociationHandle->end(); ++it) {
      reco::TrackRef itTrack  = it->val;
      edm::Ref<std::vector<Trajectory> > itraj  = it->key; // bug to find type of the key
      std::vector<TrajectoryMeasurement> TMeas=itraj->measurements();
      vector<TrajectoryMeasurement>::iterator itm;
             
      // remove tracks with high chi2
      //if(itTrack->normalizedChi2() > 5) continue;
      
      TreeTrack treetrack; // container to store tracks infos      
      treetrack.pT = itTrack->pt();
      treetrack.p = itTrack->p();
      treetrack.chi2 = itTrack->normalizedChi2();
	  treeevent->Ntracks++;
      
	  // Loop on trajectory measurements
      for (itm=TMeas.begin();itm!=TMeas.end();itm++){
	
		const TrackingRecHit* hit = &*(*itm).recHit();
		const DetId detid = hit->geographicalId();
		int subDet = detid.subdetId();
		
		treetrack.Nhits++;
		
		// For small bias scan, store only info of modules from powergroups used for the scan
		if(isSmallBiasScan_)
		if( !isUsedInSmallBiasScan(detid.rawId()) ) continue; 
		

		const TransientTrackingRecHit::ConstRecHitPointer theTTrechit = (*itm).recHit();
		const SiStripRecHit2D* SiStriphit2D = dynamic_cast<const SiStripRecHit2D*>((*theTTrechit).hit());
		const SiStripRecHit1D* SiStriphit1D = dynamic_cast<const SiStripRecHit1D*>((*theTTrechit).hit());

		const SiStripCluster* cluster = 0;

		if (SiStriphit2D!=0){
		  cluster = &*(SiStriphit2D->cluster());

		}
		if (SiStriphit1D!=0){
		  cluster = &*(SiStriphit1D->cluster());
		}


		TreeFullHit* treehit = new TreeFullHit(); // container to store hits infos
		treehit->detId = detid.rawId();

		double charge = 0;
		const TrajectoryStateOnSurface tsos = itm->updatedState();
		double angle = tsos.localDirection().theta();
                double angle_dxdz = atan2(tsos.localDirection().x(), tsos.localDirection().z());
		treehit->angle =  angle_dxdz;
		treehit->tsosx = tsos.localPosition().x();
		treehit->tsosy = tsos.localPosition().y();


		SiStripClusterInfo clusterInfo = SiStripClusterInfo( *cluster, iSetup, detid); 
		//if(cluster == 0) cout << "no cluster found " << endl;
		if(cluster!=0){

		  clusterInfo.signalOverNoise();
		  treehit->noise  = clusterInfo.noiseRescaledByGain();
		  charge = clusterInfo.charge();
		  treehit->chargeAngleCorr = charge*fabs(cos(angle));
		  treehit->width = cluster->amplitudes().size();
		  treehit->barycenter = cluster->barycenter();
		  treehit->seed = clusterInfo.maxStrip();

    	  // cluster position ----
		  DetId clusterDetId(detid);
    	  const StripGeomDetUnit * stripdet = (const StripGeomDetUnit*) theTrackerGeometry->idToDetUnit(clusterDetId);
    	  StripClusterParameterEstimator::LocalValues parameters=stripcpe.localParameters(*cluster, *stripdet);
    	  LocalPoint lp = parameters.first;
		  treehit->clusx = lp.x();

  		}

		if(abs(angle)>4) std::cout<<"ERROR : angle "<<angle<<"  evt "<<treeevent->ev_nr<<std::endl;


		//******************
		//fill hits 
		//******************
        if(charge!=0)
		{
		  if(subDet == SiStripDetId::TIB)
		  {
			// TreeFullHit inherit from TreeHit and contains more infos
			if(!fullHitInfo_TIB_) treetrack.TIB_hits.push_back( *(dynamic_cast< TreeHit* >(treehit)) );
			else treetrack.TIB_fullHits.push_back( *treehit );
		  }
		  if(subDet == SiStripDetId::TOB)
		  {
			if(!fullHitInfo_TOB_) treetrack.TOB_hits.push_back( *(dynamic_cast< TreeHit* >(treehit)) );
			else treetrack.TOB_fullHits.push_back( *treehit );
		  }
		  if(subDet == SiStripDetId::TID)
		  {
			if(!fullHitInfo_TID_) treetrack.TID_hits.push_back( *(dynamic_cast< TreeHit* >(treehit)) );
			else treetrack.TID_fullHits.push_back( *treehit );
		  }
		  if(subDet == SiStripDetId::TEC)
		  {
			if(!fullHitInfo_TEC_) treetrack.TEC_hits.push_back( *(dynamic_cast< TreeHit* >(treehit)) );
			else treetrack.TEC_fullHits.push_back( *treehit );
	      }
		}
		
		delete treehit;
		
      } // End of Loop on trajectory measurements
    
	  // store only tracks with hits stored
	  if(treetrack.TIB_hits.size() || treetrack.TIB_fullHits.size() ||
	     treetrack.TOB_hits.size() || treetrack.TOB_fullHits.size() ||
	     treetrack.TID_hits.size() || treetrack.TID_fullHits.size() ||
	     treetrack.TEC_hits.size() || treetrack.TEC_fullHits.size())
	    treeevent->tracks.push_back(treetrack);
    
      
    } // End of Loop on tracks
	
      smalltree->Fill();
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
SignalBiasScan::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SignalBiasScan::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(SignalBiasScan);
