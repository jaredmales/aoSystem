/** \file aoSystem.cpp
 * \author Jared R. Males (jaredmales@gmail.com)
 * \brief Application class and main definition for the aoSystem app.
 * \ingroup files
 *
 */

//***********************************************************************//
// Copyright 2018-2021 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of aoSystem.
//
// aoSystem is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// aoSystem is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with aoSystem.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#include <iostream>
#include <fstream>

#include <Eigen/Dense>

#include <mx/ioutils/fits/fitsFile.hpp>
#include <mx/math/constants.hpp>

#include <mx/math/func/airyPattern.hpp>

#include <mx/app/application.hpp>
#include <mx/math/ft/fftwEnvironment.hpp>

using namespace mx::math;

#include <mx/ao/analysis/aoSystem.hpp>
#include <mx/ao/analysis/aoPSDs.hpp>
#include <mx/ao/analysis/aoWFS.hpp>
#include <mx/ao/analysis/varmapToImage.hpp>
#include <mx/ao/analysis/fourierTemporalPSD.hpp>

/* git version header, generated by gengithead.sh using mxlib make system */
#include "aoSystem_git_version.h"

using namespace mx::app;

///
/**
 * Star Magnitudes:
 * - if <b>starMag</b> alone is set, then results are provided for just this one star magnitude.
 * - if <b>starMags</b>, a vector, is set, then many results are provided for each magnitude. E.g. <b>--mode</b>=ErrorBudget will produce a table.
 */
template <typename _realT>
class mxAOSystem_app : public mx::app::application
{

public:
    typedef _realT realT; ///< Real floating point type for calculations

    typedef Eigen::Array<realT, -1, -1> imageT; ///< The image type

    typedef mx::AO::analysis::aoSystem<realT, mx::AO::analysis::vonKarmanSpectrum<realT>, std::ostream> m_aosysT; ///< The AO system type.

    /// Default constructor
    mxAOSystem_app();

    /// Desctructor
    ~mxAOSystem_app();

protected:
    m_aosysT m_aosys; ///< The ao system.

    mx::AO::analysis::wfs<realT> idealWFS;                     ///< An ideal WFS
    mx::AO::analysis::pywfsUnmod<realT> unmodPyWFS;            ///< An unmodulated Pyramid WFS
    mx::AO::analysis::pywfsModAsymptotic<realT> asympModPyWFS; ///< A modulated Pyramid WFS in its asymptotic limit

    realT lam_0;

    bool m_dumpSetup{true};
    std::string m_setupOutName{"aoSystem_setup.txt"};

    std::string m_mode{"C2Var"};

    std::string wfeUnits;

    bool m_strehlOG {1};

    std::vector<realT> m_starMags;

    realT m_dfreq{0.1};
    realT m_fmax{0};
    realT k_m;
    realT k_n;
    std::string gridDir; ///< The directory for writing the grid of PSDs.
    std::string subDir;  ///< The sub-directory of gridDir where to write the analysis results.
    realT m_gfixed{0};
    int lpNc {0};                            ///< Number of linear predictor coefficients.  If <= 1 then not used.
    realT lpRegPrec {2};                 ///< The initial precision for the LP regularization algorithm.  Normal value is 2.  Higher is faster. Decrease if getting stuck in local minima.
    bool m_uncontrolledLifetimes{false}; ///< Whether or not lifetimes are calcualated for uncontrolled modes.
    int m_lifetimeTrials{0};             ///< Number of trials to use for calculating speckle lifetimes.  If 0, lifetimes are not calcualted.
    bool m_writePSDs{false};             ///< Flag controlling whether output temporal PSDs are written to disk or not.
    bool m_writeXfer{false};             ///< Flag controlling whether output transfer functions are written to disk or not.

    virtual void setupConfig();

    virtual void loadConfig();

    virtual int execute();

    int C_MapCon(const std::string &mapFile,
                 imageT &map);

    int C0Var();
    int C0Con();

    int C1Var();
    int C1Con();

    int C2Var();
    int C2Con();

    int C3Var();
    int C3Con();

    int C4Var();
    int C4Con();

    int C5Var();
    int C5Con();

    int C6Var();
    int C6Con();

    int C7Var();
    int C7Con();

    int CVarAll();

    int CConAll();

    int ErrorBudget();

    int temporalPSD();

    int temporalPSDGrid();

    int temporalPSDGridAnalyze();

    /// Output current parameters to a stream
    /** Outputs a formatted list of all current parameters.
     *
     */
    template <typename iosT>
    iosT &dumpSetup(iosT &ios /**< [in] a std::ostream-like stream. */);
};

template <typename realT>
mxAOSystem_app<realT>::mxAOSystem_app()
{
    m_configPathGlobal_env = "MXAOSYSTEM_GLOBAL_CONFIG";
    m_configPathLocal = "aoSystem.conf";

    m_requireConfigPathLocal = false;

    lam_0 = 0;

    // m_aosys.wfsBeta( idealWFS );

    wfeUnits = "rad";

    m_aosys.loadMagAOX();

    k_m = 1;
    k_n = 0;
}

template <typename realT>
mxAOSystem_app<realT>::~mxAOSystem_app()
{
}

template <typename realT>
void mxAOSystem_app<realT>::setupConfig()
{
    // App config
    config.add("mode", "m", "mode", argType::Required, "", "mode", false, "string", "Mode of calculation: ErrorBudget, C<N>Var, C<N>Con, CVarAll, CConAll,  temporalPSD, temporalPSDGrid, temporalPSDGridAnalyze");
    config.add("setupOutFile", "", "setupOutFile", argType::Required, "", "setupOutFile", false, "string", "Filename for output of setup data");
    config.add("wfeUnits", "", "wfeUnits", argType::Required, "", "wfeUnits", false, "string", "Units for WFE in ErrorBudget: rad or nm");

    // Load a model
    config.add("model", "", "model", argType::Required, "", "model", false, "string", "Model to load: Guyon2005, MagAOX, or GMagAOX");

    m_aosys.setupConfig(config);
    config.add("aosys.starMags", "", "aosys.starMags", argType::Required, "aosys", "starMags", false, "real vector", "A vector of star magnitudes");
    config.add("aosys.strehlOG", "", "aosys.strehlOG", argType::Required, "aosys", "strehlOG", false, "bool", "Flag controlling whether Strehl is used as the optical gain in error budgets. Default true.");
    // PSD Configuration
    // config.add("subTipTilt"    ,"", "subTipTilt",    argType::Required, "PSD", "subTipTilt",    false, "bool",   "If set to true, the Tip/Tilt component is subtracted from the PSD.");
    // config.add("scintillation" ,"", "scintillation", argType::Required, "PSD", "scintillation", false, "bool",   "If set to true, then scintillation is included in the PSD.");
    // config.add("component"     ,"", "component",     argType::Required, "PSD", "component",     false, "string", "Can be phase [default], amplitude, or dispersion.");

    // Temporal configuration
    config.add("fmax", "", "fmax", argType::Required, "temporal", "fmax", false, "real", "Maximum temporal frequency at which to explicitly calculate PSDs.  If 0 (default) this is based on highest wind peak.  A -17/3 power law is used above this frequency.");
    config.add("dfreq", "", "dfreq", argType::Required, "temporal", "dfreq", false, "real", "Spacing of frequencies in the analysis.");
    config.add("k_m", "", "k_m", argType::Required, "temporal", "k_m", false, "real", "The spatial frequency m index.");
    config.add("k_n", "", "k_n", argType::Required, "temporal", "k_n", false, "real", "The spatial frequency n index.");
    config.add("gridDir", "", "gridDir", argType::Required, "temporal", "gridDir", false, "string", "The directory to store the grid of PSDs.");
    config.add("subDir", "", "subDir", argType::Required, "temporal", "subDir", false, "string", "The directory to store the analysis results.");
    config.add("gfixed", "", "gfixed", argType::Required, "temporal", "gfixed", false, "float", "if > 0 then this fixed gain is used in the SI.");
    config.add("lpNc", "", "lpNc", argType::Required, "temporal", "lpNc", false, "int", "The number of linear prediction coefficients to use (if <= 1 ignored)");
    config.add("lpRegPrec", "", "lpRegPrec", argType::Required, "temporal", "lpRegPrec", false, "int", "The initial precision for the LP regularization algorithm.  Normal value is 2.  Higher is faster. Decrease if getting stuck in local minima.s to use (if <= 1 ignored)");
    config.add("uncontrolledLifetimes", "", "uncontrolledLifetimes", argType::Required, "temporal", "uncontrolledLifetimes", false, "bool", "If true, lifetimes are calculated for uncontrolled modes.  Default is false.");
    config.add("lifetimeTrials", "", "lifetimeTrials", argType::Required, "temporal", "lifetimeTrials", false, "int", "Number of trials to use for calculating speckle lifetimes.  If 0, lifetimes are not calcualted.");
    config.add("writePSDs", "", "writePSDs", argType::True, "temporal", "writePSDs", false, "bool", "Flag.  If set then output PSDs are written to disk.");
    config.add("writeXfer", "", "writeXfer", argType::True, "temporal", "writeXfer", false, "bool", "Flag.  If set then output transfer functions are written to disk.");
}

template <typename realT>
void mxAOSystem_app<realT>::loadConfig()
{
    std::vector<realT> vecTmp;

    /**********************************************************/
    /* App setup                                              */
    /**********************************************************/
    config(doHelp, "help");

    config(m_mode, "mode");

    config(wfeUnits, "wfeUnits");

    /**********************************************************/
    /* Models                                                 */
    /**********************************************************/
    // These are called before anything else, so all parameters are modifications.
    ///\todo this should actually get processed as a config file, so these settings override previous configs for proper cascading.  As it is now, if D is set in any config file, it will get overwritten, for instance.
    std::string model;
    config(model, "model");

    if(model != "")
    {
        if(model == "Guyon2005")
            m_aosys.loadGuyon2005();
        else if(model == "MagAOX")
            m_aosys.loadMagAOX();
        else if(model == "GMagAOX")
            m_aosys.loadGMagAOX();
        else
        {
            std::cerr << "Unknown model: " << model << "\n";
            doHelp = true;
            return;
        }
    }


    /**********************************************************/
    /* System                                                 */
    /**********************************************************/

    m_aosys.loadConfig(config);

    if(config.isSet("aosys.starMags"))
    {
        config(m_starMags, "aosys.starMags");
    }

    config(m_strehlOG, "aosys.strehlOG");

    /**********************************************************/
    /* Temporal PSDs                                          */
    /**********************************************************/
    config(m_fmax, "fmax");
    config(m_dfreq, "dfreq");
    config(k_m, "k_m");
    config(k_n, "k_n");

    config(gridDir, "gridDir");
    config(subDir, "subDir");
    config(m_gfixed, "gfixed");
    config(lpNc, "lpNc");
    config(lpRegPrec, "lpRegPrec");

    config(m_uncontrolledLifetimes, "uncontrolledLifetimes");
    config(m_lifetimeTrials, "lifetimeTrials");
    config(m_writePSDs, "writePSDs");
    config(m_writeXfer, "writeXfer");

    if(config.m_unusedConfigs.size() > 0)
    {
        std::cerr << "****************************************************\n";
        std::cerr << "WARNING: unrecognized config options:\n";

        for(auto it = config.m_unusedConfigs.begin(); it != config.m_unusedConfigs.end(); ++it)
        {
            std::cerr << "   " << it->second.name;
            if(config.m_sources)
                std::cerr << " [" << it->second.sources[0] << "]\n";
            else
                std::cerr << "\n";
        }

        std::cerr << "****************************************************\n";
    }

    if(config.nonOptions.size() > 0)
    {
        std::cerr << "****************************************************\n";
        std::cerr << "WARNING: unrecognized command line arguments\n";
    }
}

template <typename realT>
int mxAOSystem_app<realT>::execute()
{
    int rv;

    if(m_mode == "ErrorBudget")
    {
        rv = ErrorBudget();
    }
    else if(m_mode == "C0Var")
    {
        rv = C0Var();
    }
    else if(m_mode == "C0Con")
    {
        rv = C0Con();
    }
    else if(m_mode == "C1Var")
    {
        rv = C1Var();
    }
    else if(m_mode == "C1Con")
    {
        rv = C1Con();
    }
    else if(m_mode == "C2Var")
    {
        rv = C2Var();
    }
    else if(m_mode == "C2Con")
    {
        rv = C2Con();
    }
    else if(m_mode == "C3Var")
    {
        rv = C3Var();
    }
    else if(m_mode == "C3Con")
    {
        rv = C3Con();
    }
    else if(m_mode == "C4Var")
    {
        rv = C4Var();
    }
    else if(m_mode == "C4Con")
    {
        rv = C4Con();
    }
    else if(m_mode == "C5Var")
    {
        rv = C5Var();
    }
    else if(m_mode == "C5Con")
    {
        rv = C5Con();
    }
    else if(m_mode == "C6Var")
    {
        rv = C6Var();
    }
    else if(m_mode == "C6Con")
    {
        rv = C6Con();
    }
    else if(m_mode == "C7Var")
    {
        rv = C7Var();
    }
    else if(m_mode == "C7Con")
    {
        rv = C7Con();
    }
    else if(m_mode == "CVarAll")
    {
        rv = CVarAll();
    }
    else if(m_mode == "CConAll")
    {
        rv = CConAll();
    }
    else if(m_mode == "temporalPSD")
    {
        rv = temporalPSD();
    }
    else if(m_mode == "temporalPSDGrid")
    {
        rv = temporalPSDGrid();
    }
    else if(m_mode == "temporalPSDGridAnalyze")
    {
        rv = temporalPSDGridAnalyze();
    }
    else
    {
        std::cerr << "Unknown mode: " << m_mode << "\n\n";
        rv = -1;
    }

    if(m_dumpSetup)
    {
        std::ofstream fout;
        fout.open(m_setupOutName);

        dumpSetup(fout);

        fout.close();
    }

    if(rv < 0)
    {
        std::cerr << "An error occurred. Please use -h or --help to get help.";
        if(m_dumpSetup)
            std::cerr << "  Current settings written to " << m_setupOutName;
        std::cerr << "\n";
    }

    return rv;
}

template <typename realT>
int mxAOSystem_app<realT>::C_MapCon(const std::string &mapFile,
                                    imageT &map)
{
    imageT im, psf;

    psf.resize(map.rows(), map.cols());
    for(int i = 0; i < psf.rows(); ++i)
    {
        for(int j = 0; j < psf.cols(); ++j)
        {
            psf(i, j) = mx::math::func::airyPattern(sqrt(pow(i - floor(.5 * psf.rows()), 2) + pow(j - floor(.5 * psf.cols()), 2)));
        }
    }

    mx::AO::analysis::varmapToImage(im, map, psf);

    for(int i = 0; i < m_aosys.fit_mn_max(); ++i)
    {
        std::cout << i << " " << im(m_aosys.fit_mn_max() + 1, m_aosys.fit_mn_max() + 1 + i) << "\n";
    }

    mx::fits::fitsFile<realT> ff;
    ff.write(mapFile, im);

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::C0Var()
{
    imageT map;

    map.resize(m_aosys.fit_mn_max() * 2 + 1, m_aosys.fit_mn_max() * 2 + 1);
    m_aosys.C0Map(map, false);

    mx::fits::fitsFile<realT> ff;
    ff.write("C0Var.fits", map);

    for(int i = 0; i < m_aosys.fit_mn_max(); ++i)
    {
        std::cout << i << " " << m_aosys.C0(i, 0, false) << "\n";
    }

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::C0Con()
{
    imageT map;

    map.resize(m_aosys.fit_mn_max() * 2 + 1, m_aosys.fit_mn_max() * 2 + 1);

    m_aosys.C0Map(map, true);

    C_MapCon("C0Con.fits", map);

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::C1Var()
{
    imageT map;

    map.resize(m_aosys.fit_mn_max() * 2 + 1, m_aosys.fit_mn_max() * 2 + 1);
    m_aosys.C1Map(map, false);

    mx::fits::fitsFile<realT> ff;
    ff.write("C1Var.fits", map);

    for(int i = 0; i < m_aosys.fit_mn_max(); ++i)
    {
        std::cout << i << " " << m_aosys.C1(i, 0, false) << "\n";
    }

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::C1Con()
{
    imageT map;

    map.resize(m_aosys.fit_mn_max() * 2 + 1, m_aosys.fit_mn_max() * 2 + 1);

    m_aosys.C1Map(map, true);

    C_MapCon("C1Con.fits", map);

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::C2Var()
{
    imageT map;

    map.resize(m_aosys.fit_mn_max() * 2 + 1, m_aosys.fit_mn_max() * 2 + 1);

    m_aosys.C2Map(map, false);

    mx::fits::fitsFile<realT> ff;
    ff.write("C2Var.fits", map);

    for(int i = 0; i < m_aosys.fit_mn_max(); ++i)
    {
        std::cout << i << " " << m_aosys.C2(i, 0, false) << "\n";
    }

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::C2Con()
{
    imageT map;

    map.resize(m_aosys.fit_mn_max() * 2 + 1, m_aosys.fit_mn_max() * 2 + 1);

    m_aosys.C2Map(map, true);

    C_MapCon("C2Con.fits", map);

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::C3Var()
{
    imageT map;

    map.resize(m_aosys.fit_mn_max() * 2 + 1, m_aosys.fit_mn_max() * 2 + 1);
    m_aosys.C3Map(map, false);

    mx::fits::fitsFile<realT> ff;
    ff.write("C3Var.fits", map);

    for(int i = 0; i < m_aosys.fit_mn_max(); ++i)
    {
        std::cout << i << " " << m_aosys.C3(i, 0, false) << "\n";
    }

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::C3Con()
{
    imageT map;

    map.resize(m_aosys.fit_mn_max() * 2 + 1, m_aosys.fit_mn_max() * 2 + 1);

    m_aosys.C3Map(map, true);

    C_MapCon("C3Con.fits", map);

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::C4Var()
{
    imageT map;

    map.resize(m_aosys.fit_mn_max() * 2 + 1, m_aosys.fit_mn_max() * 2 + 1);
    m_aosys.C4Map(map, false);

    mx::fits::fitsFile<realT> ff;
    ff.write("C4Var.fits", map);

    for(int i = 0; i < m_aosys.fit_mn_max(); ++i)
    {
        std::cout << i << " " << m_aosys.C4(i, 0, false) << "\n";
    }

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::C4Con()
{
    imageT map;

    map.resize(m_aosys.fit_mn_max() * 2 + 1, m_aosys.fit_mn_max() * 2 + 1);

    m_aosys.C4Map(map, true);

    C_MapCon("C4Con.fits", map);

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::C5Var()
{
    imageT map;

    map.resize(m_aosys.fit_mn_max() * 2 + 1, m_aosys.fit_mn_max() * 2 + 1);
    m_aosys.C5Map(map, false);

    mx::fits::fitsFile<realT> ff;
    ff.write("C5Var.fits", map);

    for(int i = 0; i < m_aosys.fit_mn_max(); ++i)
    {
        std::cout << i << " " << m_aosys.C5(i, 0, false) << "\n";
    }

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::C5Con()
{
    imageT map;

    map.resize(m_aosys.fit_mn_max() * 2 + 1, m_aosys.fit_mn_max() * 2 + 1);

    m_aosys.C5Map(map, true);

    C_MapCon("C5Con.fits", map);

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::C6Var()
{
    imageT map;

    map.resize(m_aosys.fit_mn_max() * 2 + 1, m_aosys.fit_mn_max() * 2 + 1);
    m_aosys.C6Map(map, false);

    mx::fits::fitsFile<realT> ff;
    ff.write("C6Var.fits", map);

    for(int i = 0; i < m_aosys.fit_mn_max(); ++i)
    {
        std::cout << i << " " << m_aosys.C6(i, 0, false) << "\n";
    }

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::C6Con()
{
    imageT map;

    map.resize(m_aosys.fit_mn_max() * 2 + 1, m_aosys.fit_mn_max() * 2 + 1);

    m_aosys.C6Map(map, true);

    C_MapCon("C6Con.fits", map);

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::C7Var()
{
    imageT map;

    map.resize(m_aosys.fit_mn_max() * 2 + 1, m_aosys.fit_mn_max() * 2 + 1);
    m_aosys.C7Map(map, false);

    mx::fits::fitsFile<realT> ff;
    ff.write("C7Var.fits", map);

    for(int i = 0; i < m_aosys.fit_mn_max(); ++i)
    {
        std::cout << i << " " << m_aosys.C7(i, 0, false) << "\n";
    }

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::C7Con()
{
    imageT map;

    map.resize(m_aosys.fit_mn_max() * 2 + 1, m_aosys.fit_mn_max() * 2 + 1);

    m_aosys.C7Map(map, true);

    C_MapCon("C7Con.fits", map);

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::CVarAll()
{
    std::cout << "# Residual Variance Profiles\n";
    std::cout << "# sep[lam/d]    C0    C1    C2    C3    C4    C5    C6   C7\n";
    for(int i = 0; i < m_aosys.fit_mn_max(); ++i)
    {
        std::cout << i << " " << m_aosys.C0(i, 0, false) << " " << m_aosys.C1(i, 0, false) << " ";
        std::cout << m_aosys.C2(i, 0, false) << " " << m_aosys.C3(i, 0, false) << " ";
        std::cout << m_aosys.C4(i, 0, false) << " " << m_aosys.C5(i, 0, false) << " ";
        std::cout << m_aosys.C6(i, 0, false) << " " << m_aosys.C7(i, 0, false) << "\n";
    }

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::CConAll()
{
    imageT map, im0, im1, im2, im3, im4, im5, im6, im7;

    map.resize(m_aosys.fit_mn_max() * 2 + 1, m_aosys.fit_mn_max() * 2 + 1);

    imageT psf;

    psf.resize(map.rows(), map.cols());
    for(int i = 0; i < psf.rows(); ++i)
    {
        for(int j = 0; j < psf.cols(); ++j)
        {
            psf(i, j) = mx::math::func::airyPattern(sqrt(pow(i - floor(.5 * psf.rows()), 2) + pow(j - floor(.5 * psf.cols()), 2)));
        }
    }

    m_aosys.C0Map(map, true);
    mx::AO::analysis::varmapToImage(im0, map, psf);

    m_aosys.C1Map(map, true);
    mx::AO::analysis::varmapToImage(im1, map, psf);

    m_aosys.C2Map(map, true);
    mx::AO::analysis::varmapToImage(im2, map, psf);

    m_aosys.C3Map(map, true);
    mx::AO::analysis::varmapToImage(im3, map, psf);

    m_aosys.C4Map(map, true);
    mx::AO::analysis::varmapToImage(im4, map, psf);

    m_aosys.C5Map(map, true);
    mx::AO::analysis::varmapToImage(im5, map, psf);

    m_aosys.C6Map(map, true);
    mx::AO::analysis::varmapToImage(im6, map, psf);

    m_aosys.C7Map(map, true);
    mx::AO::analysis::varmapToImage(im7, map, psf);

    std::cout << "# Residual Contrast Profiles\n";
    std::cout << "# sep[lam/d]    C0    C1    C2    C3    C4    C5    C6   C7\n";
    for(int i = 0; i < m_aosys.fit_mn_max(); ++i)
    {
        std::cout << i + 1 << " " << im0(m_aosys.fit_mn_max() + 1, m_aosys.fit_mn_max() + 1 + i) << " " << im1(m_aosys.fit_mn_max() + 1, m_aosys.fit_mn_max() + 1 + i) << " ";
        std::cout << im2(m_aosys.fit_mn_max() + 1, m_aosys.fit_mn_max() + 1 + i) << " " << im3(m_aosys.fit_mn_max() + 1, m_aosys.fit_mn_max() + 1 + i) << " ";
        std::cout << im4(m_aosys.fit_mn_max() + 1, m_aosys.fit_mn_max() + 1 + i) << " " << im5(m_aosys.fit_mn_max() + 1, m_aosys.fit_mn_max() + 1 + i) << " ";
        std::cout << im6(m_aosys.fit_mn_max() + 1, m_aosys.fit_mn_max() + 1 + i) << " " << im7(m_aosys.fit_mn_max() + 1, m_aosys.fit_mn_max() + 1 + i) << "\n";
    }

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::ErrorBudget()
{
    realT units = 1;

    if(wfeUnits == "nm")
    {
        units = m_aosys.lam_sci() / (mx::math::two_pi<realT>()) / 1e-9;
    }

    if(m_starMags.size() == 0)
    {
        m_starMags.push_back(m_aosys.starMag());
    }

        std::cout << "#mag\t    d\t  bin \t n_modes \t Measurement\t   TimeDelay\t Fitting\t ChrScintOPD\t    ChrIndex\t    DispAnisoOPD    NCP    \t     Strehl\n";

        for(size_t i = 0; i < m_starMags.size(); ++i)
        {
            m_aosys.starMag(m_starMags[i]);

            if(m_strehlOG)
            {
                // Iterative optical gain
                /// \todo need upstream NCP and CP NCP and NCP NCP
                realT ncp = m_aosys.ncp_wfe(); // save ncp and then set it to zero for this part.
                m_aosys.ncp_wfe(0);

                realT lam_sci = m_aosys.lam_sci();
                m_aosys.lam_sci(m_aosys.lam_wfs());

                realT S = m_aosys.strehl();
                // std::cerr << S << "\n";

                if(S > -1) // 0.1)
                {
                    for(int s = 0; s < 4; ++s)
                    {
                        m_aosys.opticalGain(S);
                        m_aosys.optd(m_aosys.optd()); // just trigger a re-calc
                        S = m_aosys.strehl();
                        // std::cerr << S << "\n";
                    }
                }
                else
                {
                    m_aosys.opticalGain(0.1);
                }

                m_aosys.lam_sci(lam_sci);

                m_aosys.ncp_wfe(ncp);
            }

            realT S = m_aosys.strehl();
            std::cout << m_starMags[i] << "\t    ";
            std::cout << m_aosys.d_opt() << "\t  ";
            std::cout << m_aosys.bin_opt() + 1 << "\t  ";
            std::cout << m_aosys.controlledModes().sum() << "\t";
            std::cout << sqrt(m_aosys.measurementErrorTotal()) * units << "\t   ";
            std::cout << sqrt(m_aosys.timeDelayErrorTotal()) * units << "\t ";
            std::cout << sqrt(m_aosys.fittingErrorTotal()) * units << "\t ";
            std::cout << sqrt(m_aosys.chromScintOPDErrorTotal()) * units << "\t    ";
            std::cout << sqrt(m_aosys.chromIndexErrorTotal()) * units << "\t    ";
            std::cout << sqrt(m_aosys.dispAnisoOPDErrorTotal()) * units << "\t    ";
            std::cout << sqrt(m_aosys.ncpError()) * units << "\t    ";
            std::cout << S << "\n";
        }

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::temporalPSD()
{
    std::vector<realT> freq, psdOL, psdN, psdSI, psdLP;

    mx::AO::analysis::fourierTemporalPSD<realT, m_aosysT> ftPSD;
    ftPSD.m_aosys = &m_aosys;

    if(m_aosys.tauWFS() <= 0)
    {
        std::cerr << "temporalPSD: You must set tauWFS to be > 0 to specify loop frequency.\n";
        return -1;
    }

    if(m_dfreq <= 0)
    {
        std::cerr << "temporalPSD: You must set dfreq to be > 0 to specify frequency sampling.\n";
        return -1;
    }

    realT fs = 1.0 / m_aosys.tauWFS();

    mx::math::vectorScale(freq, 0.5 * fs / m_dfreq, m_dfreq, m_dfreq);
    psdOL.resize(freq.size());

    ftPSD.multiLayerPSD(psdOL, freq, k_m, k_n, 1, m_fmax);

    // Create WFS Noise PSD
    psdN.resize(freq.size());
    mx::AO::analysis::wfsNoisePSD<realT>(psdN, m_aosys.beta_p(k_m, k_n), m_aosys.Fg(), (1.0 / fs), m_aosys.npix_wfs((size_t)0), m_aosys.Fbg((size_t)0), m_aosys.ron_wfs((size_t)0));

    // optimize
    mx::AO::analysis::clAOLinearPredictor<realT> tflp;

    mx::AO::analysis::clGainOpt<realT> go_si(m_aosys.tauWFS(), m_aosys.deltaTau());
    mx::AO::analysis::clGainOpt<realT> go_lp(m_aosys.tauWFS(), m_aosys.deltaTau());

    go_si.f(freq);

    realT gmaxSI = 0;
    realT varSI;
    realT goptSI = go_si.optGainOpenLoop(varSI, psdOL, psdN, gmaxSI);

    // Only bother if number of coefficients is > 1.
    realT goptLP = -1;
    realT varLP = -1;
    if(lpNc > 1)
    {
        go_lp.f(freq);
        realT gmaxLP;
        tflp.m_precision0 = lpRegPrec;
        tflp.regularizeCoefficients(gmaxLP, goptLP, varLP, go_lp, psdOL, psdN, lpNc);
    }

    realT tauOL = 0, tauSI = 0, tauLP = 0;

    if(m_lifetimeTrials > 0)
    {
        std::vector<realT> spfreq, sppsd;
        std::vector<std::complex<realT>> ETFxn(freq.size()), NTFxn(freq.size());
        mx::sigproc::psdVarMean<mx::sigproc::psdVarMeanParams<realT>> pvm;
        mx::sigproc::psdVarMean<mx::sigproc::psdVarMeanParams<realT>> pvmLP;
        realT splifeT = 100.0;
        realT error;
        realT spvar;

        // Measure OL
        for(size_t i = 0; i < freq.size(); ++i)
        {
            ETFxn[i] = 1.0;
            NTFxn[i] = 0.0;
        }
        mx::AO::analysis::speckleAmpPSD(spfreq, sppsd, freq, psdOL, ETFxn, psdN, NTFxn, m_lifetimeTrials);
        spvar = mx::sigproc::psdVar(spfreq, sppsd);
        tauOL = pvm(error, spfreq, sppsd, splifeT) * (splifeT) / spvar;

        if(goptSI > 0)
        {
            for(size_t i = 0; i < freq.size(); ++i)
            {
                ETFxn[i] = go_si.clETF(i, goptSI);
                NTFxn[i] = go_si.clNTF(i, goptSI);
            }
        }
        else
        {
            for(size_t i = 0; i < freq.size(); ++i)
            {
                ETFxn[i] = 1;
                NTFxn[i] = 0;
            }
        }

        mx::AO::analysis::speckleAmpPSD(spfreq, sppsd, freq, psdOL, ETFxn, psdN, NTFxn, m_lifetimeTrials);
        spvar = mx::sigproc::psdVar(spfreq, sppsd);
        tauSI = pvm(error, spfreq, sppsd, splifeT) * (splifeT) / spvar;

        if(goptLP > 0)
        {
            for(size_t i = 0; i < freq.size(); ++i)
            {
                ETFxn[i] = go_lp.clETF(i, goptLP);
                NTFxn[i] = go_lp.clNTF(i, goptLP);
            }
        }
        else
        {
            for(size_t i = 0; i < freq.size(); ++i)
            {
                ETFxn[i] = 1;
                NTFxn[i] = 0;
            }
        }

        mx::AO::analysis::speckleAmpPSD(spfreq, sppsd, freq, psdOL, ETFxn, psdN, NTFxn, m_lifetimeTrials);
        spvar = mx::sigproc::psdVar(spfreq, sppsd);
        tauLP = pvmLP(error, spfreq, sppsd, splifeT) * (splifeT) / spvar;
    }

    dumpSetup(std::cout);

    std::cout << "# aoSystem single temporal PSD\n";
    std::cout << "#    var OL = " << "\n";
    if(m_lifetimeTrials > 0)
        std::cout << "#    tau OL = " << tauOL << "\n";
    std::cout << "#    opt-gain SI = " << goptSI << "\n";
    std::cout << "#    var SI = " << varSI << "\n";
    if(m_lifetimeTrials > 0)
        std::cout << "#    tau SI = " << tauSI << "\n";
    std::cout << "#    gfixed = " << m_gfixed << '\n';
    std::cout << "#    LP Num. coeff = " << lpNc << "\n";
    std::cout << "#    LP Reg. Precision = " << lpRegPrec << "\n";
    std::cout << "#    opt-gain LP = " << goptLP << "\n";
    std::cout << "#    var LP = " << varLP << "\n";
    if(m_lifetimeTrials > 0)
        std::cout << "#    tau LP = " << tauLP << "\n";

    std::cout << "####################################################################################\n";
    std::cout << "# freq       PSD-OL       PSD-N       ETF2-SI     NTF2-SI      ETF2-LP      NTF2-LP \n";

    realT ETF_SI, NTF_SI, ETF_LP = -1, NTF_LP = -1;
    for(size_t i = 0; i < freq.size(); ++i)
    {
        go_si.clTF2(ETF_SI, NTF_SI, i, goptSI);

        if(goptLP != -1)
        {
            go_lp.clTF2(ETF_LP, NTF_LP, i, goptLP);
        }

        std::cout << freq[i] << " " << psdOL[i] << " " << psdN[i] << " " << ETF_SI << " " << NTF_SI << " " << ETF_LP << " " << NTF_LP << "\n";
    }

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::temporalPSDGrid()
{
    mx::AO::analysis::fourierTemporalPSD<realT, m_aosysT> ftPSD;
    ftPSD.m_aosys = &m_aosys;

    if(gridDir == "")
    {
        std::cerr << "temporalPSDGrid: You must set gridDir.\n";
        return -1;
    }

    if(m_aosys.fit_mn_max() <= 0)
    {
        std::cerr << "temporalPSDGrid: You must set fit_mn_max to be > 0.\n";
        return -1;
    }

    if(m_dfreq <= 0)
    {
        std::cerr << "temporalPSDGrid: You must set dfreq to be > 0 to specify frequency sampling.\n";
        return -1;
    }

    if(m_aosys.tauWFS() <= 0)
    {
        std::cerr << "temporalPSDGrid: You must set tauWFS to be > 0 to specify loop frequency.\n";
        return -1;
    }

    realT fs = 1.0 / m_aosys.tauWFS();

    ftPSD.makePSDGrid(gridDir, m_aosys.fit_mn_max(), m_dfreq, fs, 0);

    return 0;
}

template <typename realT>
int mxAOSystem_app<realT>::temporalPSDGridAnalyze()
{
    mx::AO::analysis::fourierTemporalPSD<realT, m_aosysT> ftPSD;
    ftPSD.m_aosys = &m_aosys;

    if(gridDir == "")
    {
        std::cerr << "temporalPSDGridAnalyze: You must set gridDir.\n";
        return -1;
    }

    if(subDir == "")
    {
        std::cerr << "temporalPSDGridAnalyze: You must set subDir.\n";
        return -1;
    }

    if(m_aosys.fit_mn_max() <= 0)
    {
        std::cerr << "temporalPSDGridAnalyze: You must set fit_mn_max to be > 0.\n";
        return -1;
    }

    int mnCon = m_aosys.D() / m_aosys.d_min((size_t)0) / 2;

    std::vector<realT> mags;

    if(m_starMags.size() == 0)
    {
        mags = {m_aosys.starMag()};
    }
    else
    {
        mags = m_starMags;
    }

    ftPSD.m_strehlOG = m_strehlOG;

    return ftPSD.analyzePSDGrid(subDir, gridDir, m_aosys.fit_mn_max(), mnCon, m_gfixed, lpNc, lpRegPrec, mags, m_lifetimeTrials, m_uncontrolledLifetimes, m_writePSDs, m_writeXfer);
}

template <typename realT>
template <typename iosT>
iosT &mxAOSystem_app<realT>::dumpSetup(iosT &ios)
{

    ios << "# AO System Setup:\n";
    ios << "#    mode = " << m_mode << '\n';
    ios << "#    starMags = ";
    if(m_starMags.size() > 0)
    {
        for(size_t n = 0; n < m_starMags.size() - 1; ++n)
            ios << m_starMags[n] << ", ";
        ios << m_starMags.back() << '\n';
    }
    else
        ios << '\n';

    if(m_mode == "temporalPSD")
    {
        ios << "#    k_m = " << k_m << '\n';
        ios << "#    k_n = " << k_n << '\n';
    }

    if(m_mode == "temporalPSD" || m_mode == "temporalPSDGrid")
    {
        ios << "#    dfreq = " << m_dfreq << '\n';
        ios << "#    fmax = " << m_fmax << '\n';
        ios << "#    gfixed = " << m_gfixed << '\n';
        ios << "#    lpNc = " << lpNc << '\n';
        ios << "#    lpRegPrec = " << lpRegPrec << '\n';
    }

    m_aosys.dumpAOSystem(ios);

    ios << "#    aoSystem branch = " << AOSYSTEM_GIT_BRANCH << '\n';
    ios << "#    aoSystem sha1 = " << AOSYSTEM_GIT_CURRENT_SHA1;
    if(AOSYSTEM_GIT_REPO_MODIFIED)
        ios << " (modified)";
    ios << '\n';

    return ios;
}

int main(int argc, char **argv)
{
    mx::math::ft::fftwEnvironment<double, false> fftwEnv;

    mxAOSystem_app<double> m_aosysA;

    return m_aosysA.main(argc, argv);
}
