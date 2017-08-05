


#include <iostream>
#include <fstream>

#include <Eigen/Dense>


#include <mx/improc/fitsFile.hpp>
#include <mx/math/func/airyPattern.hpp>

#define MX_APP_DEFAULT_configPathGlobal_env "MXAOSYSTEM_GLOBAL_CONFIG"
#define MX_APP_DEFAULT_configPathLocal "aoSystem.conf"

#include <mx/app/application.hpp>
#include <mx/ao/analysis/aoSystem.hpp>
#include <mx/ao/analysis/aoPSDs.hpp>
#include <mx/ao/analysis/aoWFS.hpp>
#include <mx/ao/analysis/varmapToImage.hpp>
#include <mx/ao/analysis/fourierTemporalPSD.hpp>
///
/**
  * Star Magnitudes:
  * - if <b>starMag</b> alone is set, then results are provided for just this one star magnitude.
  * - if <b>starMags</b>, a vector, is set, then many results are provided for each magnitude. E.g. <b>--mode</b>=ErrorBudget will produce a table.
  */  
template<typename _realT>
class mxAOSystem_app : public mx::application
{

public:
   typedef _realT realT;

   typedef Eigen::Array<realT, -1,-1> imageT;
   
   typedef mx::AO::aoSystem<realT, mx::AO::vonKarmanSpectrum<realT>> aosysT;

   mxAOSystem_app();
   
   ~mxAOSystem_app();
   
protected:
   
   aosysT aosys;
   
   mx::AO::beta_p::wfs<realT> idealWFS;
   mx::AO::beta_p::pywfsUnmod<realT> unmodPyWFS;
   mx::AO::beta_p::pywfsModAsymptotic<realT> asympModPyWFS;
   
   
   realT lam_0;
   
   bool dumpSetup;
   std::string setupOutName;
   
   std::string mode;
   
   std::string wfeUnits;
   
   int mnMap;
   
   std::vector<realT> starMags;
   
   realT dfreq;
   realT kmax;
   realT k_m;
   realT k_n;
   std::string gridDir; ///<The directory for writing the grid of PSDs.
   
   virtual void setupConfig();

   virtual void loadConfig();
   
   virtual int execute();
   
   int C_MapCon( const std::string & mapFile,
                 imageT & map 
               );
              
   int C0Raw();
   int C0Map();
   
   int C1Raw();
   int C1Map();
   
   int C2Raw();
   int C2Map();
   
   int C4Raw();
   int C4Map();
   
   int C6Raw();
   int C6Map();
   
   int C7Raw();
   int C7Map();
   
   int CAllRaw();
   
   int ErrorBudget();
   
   int Strehl();
   
   int temporalPSD();
   
   int temporalPSDGrid();
};

template<typename realT>
mxAOSystem_app<realT>::mxAOSystem_app()
{
   lam_0 = 0;
   
   dumpSetup = true;
   setupOutName = "mxAOAnalysisSetup.txt";
   
   mode = "C2Raw";
   
   aosys.wfsBeta( idealWFS );
   
   wfeUnits = "rad";
   
   aosys.loadMagAOX(); 
   
   mnMap = 50;
   
   dfreq = 0.1;
   kmax = 0;
   k_m = 1;
   k_n = 0;

}

template<typename realT>
mxAOSystem_app<realT>::~mxAOSystem_app()
{
}   

template<typename realT>
void mxAOSystem_app<realT>::setupConfig()
{
   //App config
   config.add("mode"        ,"m", "mode" , mx::argType::Required, "", "mode",     false,  "string", "Mode of calculation: C2Raw, C2Map, ErrorBudget, Strehl");
   config.add("setupOutFile"        ,"", "setupOutFile" , mx::argType::Required, "", "setupOutFile", false, "string", "Filename for output of setup data");

   config.add("wfeUnits"        ,"", "wfeUnits" , mx::argType::Required, "", "wfeUnits", false, "string", "Units for WFE in ErrorBudget: rad or nm");

   config.add("mnMap"        ,"", "mnMap" , mx::argType::Required, "", "mnMap",     false,  "string", "Maximum spatial frequency index to include in maps.");
   
   //Load a model
   config.add("model"        ,"", "model" , mx::argType::Required, "", "model", false, "string", "Model to load: Guyon2005, MagAOX, or GMagAOX");
   
   //Atmosphere configuration
   config.add("lam_0"        ,"", "lam_0" , mx::argType::Required, "atmosphere", "lam_0",        false, "real", "The reference wavlength for r_0 [m]");
   config.add("r_0"          ,"", "r_0"   , mx::argType::Required, "atmosphere", "r_0",          false, "real", "Fried's parameter [m]");
   config.add("L_0"          ,"", "L_0"   , mx::argType::Required, "atmosphere", "L_0",          false, "real", "Outer scale [m]");
   config.add("layer_Cn2"    ,"", ""      , mx::argType::None,     "atmosphere", "layer_Cn2",    false, "real vector", "Layer Cn^2");  
   config.add("layer_v_wind" ,"", ""      , mx::argType::None,     "atmosphere", "layer_v_wind", false, "real vector", "Layer wind speeds [m/s]");
   config.add("layer_dir"    ,"", ""      , mx::argType::None,     "atmosphere", "layer_dir",    false, "real vector", "Layer wind directions [rad]");
   config.add("layer_z"      ,"", ""      , mx::argType::None,     "atmosphere", "layer_z",      false, "real vector", "layer heights [m]");
   config.add("h_obs"        ,"", ""      , mx::argType::None,     "atmosphere", "h_obs",        false, "real", "height of observatory [m]");
   config.add("H"            ,"", ""      , mx::argType::None,     "atmosphere", "H",            false, "real", "atmospheric scale heights [m]");
   config.add("v_wind"       ,"", "v_wind", mx::argType::Required, "atmosphere", "v_wind",       false, "real", "Mean windspeed (5/3 momement), rescales layers [m/s]");
   config.add("z_mean"       ,"", "z_mean", mx::argType::Required, "atmosphere", "z_mean",       false, "real", "Mean layer height (5/3 momemnt), rescales layers [m/s]");
   
   //AO System configuration
   config.add("wfs"          ,"", "wfs"        , mx::argType::Required, "system", "wfs",     false, "string", "The WFS type: idealWFS, unmodPyWFS, asympModPyWFS");
   config.add("D"            ,"", "D"          , mx::argType::Required, "system", "D",           false, "real", "The telescope diameter [m]");
   config.add("d_min"        ,"", "d_min"      , mx::argType::Required, "system", "d_min",       false, "real", "The minimum actuator spacing [m]");
   config.add("F0"           ,"", "F0"         , mx::argType::Required, "system", "F0",          false, "real", "Zero-mag photon flux, [photons/sec]");     
   config.add("lam_wfs"      ,"", "lam_wfs"    , mx::argType::Required, "system", "lam_wfs",     false, "real", "WFS wavelength [m]" );
   config.add("npix_wfs"     ,"", "npix_wfs"   , mx::argType::Required, "system", "npix_wfs",    false, "real", "The number of pixels in the WFS");
   config.add("ron_wfs"      ,"", "ron_wfs"    , mx::argType::Required, "system", "ron_wfs",     false, "real", "WFS readout noise [photons/read]");
   config.add("Fbg"          ,"", "Fbg"        , mx::argType::Required, "system", "Fbg",         false, "real", "Background counts, [counts/pix/sec]");
   config.add("minTauWFS"    ,"", "minTauWFS"  , mx::argType::Required, "system", "minTauWFS",   false, "real", "Minimum WFS integration time [s]");
   config.add("deltaTau"     ,"", "deltaTau"   , mx::argType::Required, "system", "deltaTau",    false, "real", "Loop delay [s]");
   config.add("lam_sci"      ,"", "lam_sci"    , mx::argType::Required, "system", "lam_sci",     false, "real", "Science wavelength [m]");
   config.add("zeta"         ,"", "zeta"       , mx::argType::Required, "system", "zeta",        false, "real", "Zenith distance [rad]");
   config.add("fit_mn_max"   ,"", "fit_mn_max" , mx::argType::Required, "system", "fit_mn_max",  false, "real", "Maximum spatial frequency index to use for analysis");
   config.add("ncp_wfe"      ,"", "ncp_wfe"    , mx::argType::Required, "system", "ncp_wfe",     false, "real", "NCP WFE between 1 lambda/D and fit_mn_max [rad^2]");
   config.add("ncp_alpha"    ,"", "ncp_alpha"  , mx::argType::Required, "system", "ncp_alpha",   false, "real", "PSD index for NCP WFE");
   config.add("starMag"      ,"", "starMag"    , mx::argType::Required,  "system", "starMag",     false, "real", "Star magnitude");
   config.add("starMags"     ,"", "starMags"    , mx::argType::Required,  "system", "starMags",     false, "real", "A vector of star magnitudes");
   
   //Temporal configuration
   config.add("kmax"     ,"", "kmax"    , mx::argType::Required,  "temporal", "kmax",     false, "real", "Maximum frequency at which to explicitly calculate PSDs.");
   config.add("dfreq"     ,"", "dfreq"    , mx::argType::Required,  "temporal", "dfreq",     false, "real", "Spacing of frequencies in the analysis.");
   config.add("k_m"     ,"", "k_m"    , mx::argType::Required,  "temporal", "k_m",     false, "real", "The spatial frequency m index.");
   config.add("k_n"     ,"", "k_n"    , mx::argType::Required,  "temporal", "k_n",     false, "real", "The spatial frequency n index.");
   config.add("gridDir"     ,"", "gridDir"    , mx::argType::Required,  "temporal", "gridDir",     false, "string", "The directory to store the grid of PSDs.");
   
   
   
}

template<typename realT>
void mxAOSystem_app<realT>::loadConfig()
{
   realT tmp;
   std::vector<realT> vecTmp;
   
   
   /**********************************************************/
   /* App setup                                              */
   /**********************************************************/
   config(doHelp, "help");
      
   config(mode, "mode");
   
   config(wfeUnits, "wfeUnits");
   
   config(mnMap, "mnMap");
   
   /**********************************************************/
   /* Models                                                 */
   /**********************************************************/
   //These are called before anything else, so all parameters are modifications.
   std::string model;
   config(model, "model");
   
   if(model != "")
   {
      if( model == "Guyon2005" ) aosys.loadGuyon2005(); 
      else if( model == "MagAOX" ) aosys.loadMagAOX();
      else if( model == "GMagAOX" ) aosys.loadGMagAOX();
      else
      {
         std::cerr << "Unknown model: " << model << "\n";
         return;
      }
   }
   
   /**********************************************************/
   /* Atmosphere                                             */
   /**********************************************************/
   
   //The order of l0, Cn2, and r_0 is so that r_0 overrides the value set with Cn2 if l0 != 0.
   //lam_0 comes first because it calibrates r0 and Cn2
   config(lam_0, "lam_0");
   
   //layer Cn2
   if( config.isSet("layer_Cn2") )
   {
      aosys.atm.layer_Cn2(config.get<std::vector<realT>>("layer_Cn2"), lam_0);
   }
   
   //r_0
   if(config.isSet("r_0") )
   {
      aosys.atm.r_0(config.get<realT>("r_0"), lam_0);
   }
   
   //Outer scale
   if(config.isSet("L_0") )
   {
      aosys.atm.L_0(config.get<realT>("L_0"));
   }
   
   //layer winds
   if( config.isSet("layer_v_wind") )
   {
      aosys.atm.layer_v_wind(config.get<std::vector<realT>>("layer_v_wind"));
   }
   
   //layer direction
   if( config.isSet("layer_dir") )
   {
      aosys.atm.layer_dir(config.get<std::vector<realT>>("layer_dir"));
   }
 
   if( config.isSet("layer_z") )
   {
      aosys.atm.layer_z(config.get<std::vector<realT>>("layer_z"));
   }
   

   //v_wind --> this rescales layer_v_wind
   if(config.isSet("v_wind") )
   {
      aosys.atm.v_wind(config.get<realT>("v_wind"));
   }

   //z_mean --> this rescales layer_z
   if(config.isSet("z_mean") )
   {
      aosys.atm.z_mean(config.get<realT>("z_mean"));
   }

   /**********************************************************/
   /* System                                                 */
   /**********************************************************/
   //WFS 
   if(config.isSet("wfs"))
   {
      std::string wfs;
      config(wfs, "wfs");
      
      if(wfs == "ideal")
      {
         aosys.wfsBeta(idealWFS);// = &idealWFS;
      }
      else if(wfs == "unmodPyWFS")
      {
         aosys.wfsBeta(unmodPyWFS);// = &unmodPyWFS;
      }
      else if(wfs == "asympModPyWFS")
      {
         aosys.wfsBeta(asympModPyWFS); // = &asympModPyWFS;
      }
      else
      {
         std::cerr << "Unkown WFS type\n";
         exit(-1);
      }
   }
   
   //diameter
   if(config.isSet("D") )
   {
      aosys.D(config.get<realT>("D"));
   }
   
   //d_min
   if(config.isSet("d_min") )
   {
      aosys.d_min(config.get<realT>("d_min"));
   }
      
   //F0
   if(config.isSet("F0") )
   {
      aosys.F0(config.get<realT>("F0"));
   }
   
   
   //lam_wfs
   if(config.isSet("lam_wfs") )
   {
      aosys.lam_wfs(config.get<realT>("lam_wfs"));
   }
   
   //npix_wfs
   if(config.isSet("npix_wfs") )
   {
      aosys.npix_wfs(config.get<realT>("npix_wfs"));
   }
      
   //ron_wfs
   if(config.isSet("ron_wfs") )
   {
      aosys.ron_wfs(config.get<realT>("ron_wfs"));
   }
      
   //Fbg
   if(config.isSet("Fbg") )
   {
      aosys.Fbg(config.get<realT>("Fbg"));
   }
   
   
   //minTauWFS
   if(config.isSet("minTauWFS") )
   {
      aosys.minTauWFS(config.get<realT>("minTauWFS"));
   }
      
   //deltaTau
   if(config.isSet("deltaTau") )
   {
      aosys.deltaTau(config.get<realT>("deltaTau"));
   }
      
   //lam_sci
   if(config.isSet("lam_sci") )
   {
      aosys.lam_sci(config.get<realT>("lam_sci"));
   }
   
   //zeta
   if(config.isSet("zeta") )
   {
      aosys.zeta(config.get<realT>("zeta"));
   }
   
   //fit_mn_max
   if(config.isSet("fit_mn_max") )
   {
      aosys.fit_mn_max(config.get<realT>("fit_mn_max"));
   }
   
   //ncp_wfe
   if(config.isSet("ncp_wfe") )
   {
      aosys.ncp_wfe(config.get<realT>("ncp_wfe"));
   }
    
   //ncp_alpha
   if(config.isSet("ncp_alpha") )
   {
      aosys.ncp_alpha(config.get<realT>("ncp_alpha"));
   }
      
   //star_mag
   if(config.isSet("starMag") )
   {
      aosys.starMag(config.get<realT>("starMag"));
   }
   
   if( config.isSet("starMags") )
   {
      starMags = config.get<std::vector<realT>>("starMags");
   }
   
   /**********************************************************/
   /* Temporal PSDs                                          */
   /**********************************************************/
   config.get(kmax, "kmax");
   config.get(dfreq, "dfreq");
   config.get(k_m, "k_m");
   config.get(k_n, "k_n");
   
   config.get(gridDir, "gridDir");
   
}





template<typename realT>
int mxAOSystem_app<realT>::execute()
{
   int rv;
   
   if(mode == "C0Raw")
   {
      rv = C0Raw();
   }
   else if (mode == "C0Map")
   {
      rv = C0Map();
   }
   else if(mode == "C1Raw")
   {
      rv = C1Raw();
   }
   else if (mode == "C1Map")
   {
      rv = C1Map();
   }
   else if(mode == "C2Raw")
   {
      rv = C2Raw();
   }
   else if (mode == "C2Map")
   {
      rv = C2Map();
   }
   else if(mode == "C4Raw")
   {
      rv = C4Raw();
   }
   else if (mode == "C4Map")
   {
      rv = C4Map();
   }
   else if(mode == "C6Raw")
   {
      rv = C6Raw();
   }
   else if (mode == "C6Map")
   {
      rv = C6Map();
   }
   else if(mode == "C7Raw")
   {
      rv = C7Raw();
   }
   else if (mode == "C7Map")
   {
      rv = C7Map();
   }
   else if (mode == "CAllRaw")
   {
      rv = CAllRaw();
   }
   else if (mode == "ErrorBudget")
   {
      rv = ErrorBudget();
   }
   else if (mode == "Strehl")
   {
      rv = Strehl();
   }
   else if (mode == "temporalPSD")
   {
      rv = temporalPSD();
   }
   else if (mode == "temporalPSDGrid")
   {
      rv = temporalPSDGrid();
   }
   else
   {
      std::cerr << "Unknown mode: " << mode << "\n";
      rv = -1;
   }
   
   
   
   if(dumpSetup && rv == 0)
   {
      std::ofstream fout;
      fout.open(setupOutName);
      
      aosys.dumpAOSystem(fout);
      
      fout.close();
   }
   
   return rv;
}

template<typename realT>
int mxAOSystem_app<realT>::C_MapCon( const std::string & mapFile,
                                     imageT & map
                                   )
{
   imageT im, psf;
   
   psf.resize(map.rows(),map.cols());
   for(int i=0;i<psf.rows();++i)
   {
      for(int j=0;j<psf.cols();++j)
      {
         psf(i,j) = mx::math::func::airyPattern(sqrt( pow( i-floor(.5*psf.rows()),2) + pow(j-floor(.5*psf.cols()),2)));
      }
   }
   
   mx::AO::varmapToImage(im, map, psf);
   
   for(int i=0; i< mnMap; ++i)
   {
      std::cout << i << " " << im( mnMap+1, mnMap+1 + i) << "\n";
   }
   
   mx::improc::fitsFile<realT> ff;
   ff.write(mapFile, im);

   
}

template<typename realT>
int mxAOSystem_app<realT>::C0Raw()
{
   for(int i=0;i< aosys.fit_mn_max(); ++i)
   {
      std::cout << i << " " << aosys.C0(i,0, false) << "\n";
   }
}

template<typename realT>
int mxAOSystem_app<realT>::C0Map()
{
   imageT map;
   
   map.resize( mnMap*2+1, mnMap*2 + 1);
   
   aosys.C0Map(map);
   
   C_MapCon("C0Map.fits", map);
   
}

template<typename realT>
int mxAOSystem_app<realT>::C1Raw()
{
   for(int i=0;i< aosys.fit_mn_max(); ++i)
   {
      std::cout << i << " " << aosys.C1(i,0, false) << "\n";
   }
}

template<typename realT>
int mxAOSystem_app<realT>::C1Map()
{
   imageT map;
   
   map.resize( mnMap*2+1, mnMap*2 + 1);
   
   aosys.C1Map(map);
   
   C_MapCon("C1Map.fits", map);
   
}

template<typename realT>
int mxAOSystem_app<realT>::C2Raw()
{
   for(int i=0;i< aosys.fit_mn_max(); ++i)
   {
      std::cout << i << " " << aosys.C2(i,0, false) << "\n";
   }
}

template<typename realT>
int mxAOSystem_app<realT>::C2Map()
{
   imageT map;
   
   map.resize( mnMap*2+1, mnMap*2 + 1);
   
   aosys.C2Map(map);
   
   C_MapCon("C2Map.fits", map);
   
}

template<typename realT>
int mxAOSystem_app<realT>::C4Raw()
{
   for(int i=0;i< aosys.fit_mn_max(); ++i)
   {
      std::cout << i << " " << aosys.C4(i,0, false) << "\n";
   }
}

template<typename realT>
int mxAOSystem_app<realT>::C4Map()
{
   imageT map;
   
   map.resize( mnMap*2+1, mnMap*2 + 1);
   
   aosys.C4Map(map);
   
   C_MapCon("C4Map.fits", map);
   
}

template<typename realT>
int mxAOSystem_app<realT>::C6Raw()
{
   for(int i=0;i< aosys.fit_mn_max(); ++i)
   {
      std::cout << i << " " << aosys.C6(i,0, false) << "\n";
   }
}

template<typename realT>
int mxAOSystem_app<realT>::C6Map()
{
   imageT map;
   
   map.resize( mnMap*2+1, mnMap*2 + 1);
   
   aosys.C6Map(map);
   
   C_MapCon("C6Map.fits", map);
   
}


template<typename realT>
int mxAOSystem_app<realT>::C7Raw()
{
   for(int i=0;i< aosys.fit_mn_max(); ++i)
   {
      std::cout << i << " " << aosys.C7(i,0, false) << "\n";
   }
}

template<typename realT>
int mxAOSystem_app<realT>::C7Map()
{
   imageT map;
   
   map.resize( mnMap*2+1, mnMap*2 + 1);
   
   aosys.C7Map(map);
   
   C_MapCon("C7Map.fits", map);
   
}

template<typename realT>
int mxAOSystem_app<realT>::CAllRaw()
{
   for(int i=0;i< aosys.fit_mn_max(); ++i)
   {
      std::cout << i << " " << aosys.C0(i,0, false) << " " << aosys.C1(i,0, false) << " " << aosys.C2(i,0, false) << " " << aosys.C4(i,0, false);
      std::cout << " " << aosys.C6(i,0, false) << " " << aosys.C7(i,0, false) << "\n";
   }
}


template<typename realT>
int mxAOSystem_app<realT>::ErrorBudget()
{
   realT units = 1;
   
   if(wfeUnits == "nm")
   {
      std::cerr << aosys.lam_sci() << "\n";
      units = aosys.lam_sci() / (2.0*pi<realT>()) / 1e-9;
   }
   
   if(starMags.size() == 0)
   {
      std::cout << "Measurement: " << sqrt(aosys.measurementError())*units << "\n";
      std::cout << "Time-delay:  " << sqrt(aosys.timeDelayError())*units << "\n";
      std::cout << "Fitting:     " << sqrt(aosys.fittingError())*units << "\n";
      std::cout << "NCP error:   " << sqrt(aosys.ncpError())*units << "\n";
      std::cout << "Strehl:      " << aosys.strehl() << "\n";
   }
   else
   {
      std::cout << "#mag     Measurement     Time-delay      Fitting    Chr-Scint-OPD      Chr-Index   Disp-Ansio-OPD  NCP-error         Strehl\n";
      
      for(int i=0; i< starMags.size(); ++i)
      {
         aosys.starMag(starMags[i]);
         std::cout << starMags[i] << "\t    ";
         std::cout << sqrt(aosys.measurementError())*units << "\t   ";
         std::cout << sqrt(aosys.timeDelayError())*units << "\t ";
         std::cout << sqrt(aosys.fittingError())*units << "\t ";
         std::cout << sqrt(aosys.chromScintOPDError())*units << "\t    ";
         std::cout << sqrt(aosys.chromIndexError())*units << "\t\t    ";
         std::cout << sqrt(aosys.dispAnisoOPDError())*units << "\t    ";
         std::cout << sqrt(aosys.ncpError())*units << "\t\t";
         std::cout << aosys.strehl() << "\n";
      }
   }
      
      
   return 0;
}

template<typename realT>
int mxAOSystem_app<realT>::Strehl()
{
   std::cout << aosys.strehl() << "\n";
   
   return 0;
}

template<typename realT>
int mxAOSystem_app<realT>::temporalPSD()
{
   std::vector<realT> freq, psd;

   mx::AO::fourierTemporalPSD<realT, aosysT> ftPSD;
   ftPSD._aosys = &aosys;
   
   if(aosys.minTauWFS() <= 0)
   {
      std::cerr << "temporalPSD: You must set minTauWFS to be > 0 to specify loop frequency.\n";
      return -1;
   }
   
   if(dfreq <= 0)
   {
      std::cerr << "temporalPSD: You must set dfreq to be > 0 to specify frequency sampling.\n";
      return -1;
   }
   
   realT fs = 1.0/aosys.minTauWFS();
   
   mx::vectorScale(freq, 0.5*fs/dfreq, dfreq, dfreq);
   psd.resize(freq.size());
   
   ftPSD.multiLayerPSD( psd, freq, k_m, k_n, 1, kmax);
   
   
   for(int i=0; i < freq.size(); ++i)
   {
      std::cout << freq[i] << " " << psd[i] << "\n";
   }
   
   return 0;
}

template<typename realT>
int mxAOSystem_app<realT>::temporalPSDGrid()
{
   mx::AO::fourierTemporalPSD<realT, aosysT> ftPSD;
   ftPSD._aosys = &aosys;
   
   if(gridDir == "")
   {
      std::cerr << "temporalPSDGrid: You must set gridDir.\n";
      return -1;
   }
   
   if(aosys.fit_mn_max() <= 0)
   {
      std::cerr << "temporalPSDGrid: You must set fit_mn_max to be > 0.\n";
      return -1;
   }
   

   if(dfreq <= 0)
   {
      std::cerr << "temporalPSDGrid: You must set dfreq to be > 0 to specify frequency sampling.\n";
      return -1;
   }
   
   if(aosys.minTauWFS() <= 0)
   {
      std::cerr << "temporalPSDGrid: You must set minTauWFS to be > 0 to specify loop frequency.\n";
      return -1;
   }
   
   realT fs = 1.0/aosys.minTauWFS();
   
   ftPSD.makePSDGrid( gridDir, aosys.fit_mn_max(), dfreq, fs, 0);
   
   return 0;
}


int main(int argc, char ** argv)
{
 
   mxAOSystem_app<double> aosysA;
   
   aosysA.main(argc, argv);
   

   
}

