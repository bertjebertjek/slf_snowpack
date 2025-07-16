/*
 *  SNOWPACK stand-alone
 *
 *  Copyright WSL Institute for Snow and Avalanche Research SLF, DAVOS, SWITZERLAND
*/
/*  This file is part of Snowpack.
    Snowpack is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snowpack is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snowpack.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <snowpack/SnowDrift.h>
#include <snowpack/Hazard.h>
#include <snowpack/Utils.h>

#include <vector>
#include <assert.h>

using namespace mio;
using namespace std;

/************************************************************
 * static section                                           *
 ************************************************************/

///Deviation from geometrical factors defined by Schmidt
const double SnowDrift::schmidt_drift_fudge = 1.0;

///Enables erosion notification
const bool SnowDrift::msg_erosion = false;


/************************************************************
 * non-static section                                       *
 ************************************************************/



static bool get_bool(const SnowpackConfig& cfg, const std::string& key, const std::string& section)
{
	bool value;
	cfg.getValue(key, section, value);
	return value;
}

static double get_sn_dt(const SnowpackConfig& cfg) 
{
	//Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
	const double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Snowpack");
	return M_TO_S(calculation_step_length);
}

SnowDrift::SnowDrift(const SnowpackConfig& cfg) : saltation(cfg),
                     enforce_measured_snow_heights( get_bool(cfg, "ENFORCE_MEASURED_SNOW_HEIGHTS", "Snowpack") ), snow_redistribution( get_redistribution(cfg) ), 
					 snow_erosion( get_erosion(cfg) ), alpine3d( get_bool(cfg, "ALPINE3D", "SnowpackAdvanced") ),
                     sn_dt( get_sn_dt(cfg)), erosion_limit( get_erosion_limit(cfg) )  {}

/**
 * @brief Computes the local mass flux of snow
 * @bug Contribution from suspension not considered yet!
 * @param *Edata - Single Element
 * @param ustar Shear wind velocity (m s-1)
 * @param angle Slope angle (deg)
 * @return Saltation mass flux (kg m-1 s-1)
 */
double SnowDrift::compMassFlux(const ElementData& Edata, const double& ustar, const double& slope_angle) const
{
	// Compute basic quantities that are needed: friction velocity, z0, threshold vw
	// For now assume logarithmic wind profile; TODO change this later
	const double weight = 0.02 * Constants::density_ice * (Edata.sp + 1.) * Constants::g * MM_TO_M(Edata.rg);
	// weight = Edata.Rho*(Edata.sp + 1.)*Constants::g*MM_TO_M(Edata.rg);
	const double sig = 300.;
	const double binding = 0.0015 * sig * Edata.N3 * Optim::pow2(Edata.rb/Edata.rg);
	const double tau_thresh = SnowDrift::schmidt_drift_fudge * (weight + binding);  // Original value for fudge: 1. (Schmidt)
	//const double ustar_thresh = sqrt(tau_thresh / Constants::density_air);
	const double tau = Constants::density_air * Optim::pow2(ustar);

	// First, look whether there is any transport at all: use formulation of Schmidt
	if ( tau_thresh > tau ) {
		return (0.0);
	}

	// Compute the saltation mass flux (Doorschot or Sorenson)
	double Qsalt = 0., Qsusp = 0., c_salt; // The mass fluxes in saltation and suspension (kg m-1 s-1)
	if (!saltation.compSaltation(tau, tau_thresh, slope_angle, MM_TO_M(2.*Edata.rg), Qsalt, c_salt)) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "Saltation computation failed");
		throw IOException("Saltation computation failed", AT);
	}

	if (Qsalt > 0.) {
		Qsusp = 0.; // TODO What about sd_IntegrateFlux(ustar, ustar_thresh, Qsalt, z0); ???
	} else {
		Qsalt = 0.;
		Qsusp = 0.;
	}

	return (Qsalt + Qsusp);
}


/**
 * @brief Erodes Elements from the top and computes the associated mass flux
 * Even so the code is quite obscure, it should cover all of the following cases:
 * -# Externally provided eroded mass (for example, by Alpine3D); see parameter forced_massErode
 * -# SNOW_REDISTRIBUTION is true: using vw_drift to erode the snow surface on the windward virtual slope
 * -# SNOW_EROSION is true: using vw to erode the snow surface at the main station (flat field or slope).
 * -# SNOW_EROSION is true and SNOW_REDISTRIBUTION is false: using vw_drift (if available and larger than vw)
             to do virtual erosion either on flat field (#slopes == 1) or on windward slope (#slopes > 1).
 * 	@note However, erosion will also be controlled by mH and thus a measured snow depth (HS1) is required
 * 	@note If either measured snow depth is missing or the conditions for a real erosion are not fulfilled,
 *          the possibility of a virtual erosion will be considered using the ErosionLevel marker (virtual erodible layer).
 * @param Mdata
 * @param Xdata
 * @param Sdata
 * @param forced_massErode if greater than 0, force the eroded mass to the given value (instead of computing it)
*/
void SnowDrift::compSnowDrift_deprecated(const CurrentMeteo& Mdata, SnowStation& Xdata, SurfaceFluxes& Sdata, double& forced_massErode) const
{
	size_t nE = Xdata.getNumberOfElements();
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;

	const bool no_snow = ((nE < Xdata.SoilNode+1) || (EMS[nE-1].theta[SOIL] > 0.));
	const bool no_wind_data = (Mdata.vw_drift == mio::IOUtils::nodata);
	if (no_snow || no_wind_data) {
		Xdata.ErosionMass = 0.;
		if (no_snow) {
			Xdata.ErosionLevel = Xdata.SoilNode;
			Sdata.drift = 0.;
		} else
			Sdata.drift = Constants::undefined;
		return;
	}

	// Real erosion either on windward virtual slope, from Alpine3D, or at main station.
	// At main station, measured snow depth controls whether erosion is possible or not if measured snow depth is provided
	const bool windward = !alpine3d && (snow_redistribution == "TRUE") && Xdata.windward; // check for windward virtual slope
	const bool erosion = (snow_erosion == "TRUE") && (Xdata.mH > (Xdata.Ground + Constants::eps)) && ((Xdata.mH + 0.02) < Xdata.cH);
	const double ustar_max = (Mdata.vw>0.1) ? Mdata.ustar * Mdata.vw_drift / Mdata.vw : 0.; // Scale Mdata.ustar

	if (windward || alpine3d || erosion) {
		double massErode=0.; // Mass loss due to erosion
		if (fabs(forced_massErode) > Constants::eps2) {
			massErode = std::max(0., -forced_massErode); //negative mass is erosion
		} else {
			try {
				if (enforce_measured_snow_heights && !windward)
					Sdata.drift = compMassFlux(EMS[nE-1], Mdata.ustar, Xdata.meta.getSlopeAngle()); // kg m-1 s-1, erosion at main station, local vw && nE-1
				else
					Sdata.drift = compMassFlux(EMS[nE-1], ustar_max, Xdata.meta.getSlopeAngle()); // kg m-1 s-1, windward slope && vw_drift && nE-1
			} catch(const exception&) {
					prn_msg(__FILE__, __LINE__, "err", Mdata.date, "Cannot compute mass flux of drifting snow!");
					throw;
			}
			massErode = Sdata.drift * sn_dt / Hazard::typical_slope_length; // Convert to eroded snow mass in kg m-2
		}
		unsigned int nErode=0; // number of eroded elements
		if (massErode >= 0.95 * EMS[nE-1].M) {
			// Erode at most one element with a maximal error of +- 5 % on mass ...
			if (windward)
				Xdata.rho_hn = EMS[nE-1].Rho;
			nE--;
			Xdata.cH -= EMS[nE].L;
			NDS[nE].hoar = 0.;
			Xdata.ErosionMass = EMS[nE].M;
			Xdata.ErosionLevel = std::min(nE-1, Xdata.ErosionLevel);
			nErode++;
			massErode -= EMS[nE].M;
			forced_massErode = -massErode;
		} else if (massErode > Constants::eps) { // ... or take away massErode from top element - partial real erosion
			if (fabs(EMS[nE-1].L * EMS[nE-1].Rho - EMS[nE-1].M) > 0.001) {
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date, "Inconsistent Mass:%lf   L*Rho:%lf", EMS[nE-1].M,EMS[nE-1].L*EMS[nE-1].Rho);
				EMS[nE-1].M = EMS[nE-1].L * EMS[nE-1].Rho;
				assert(EMS[nE-1].M>=0.); //mass must be positive
			}
			if (windward)
				Xdata.rho_hn = EMS[nE-1].Rho; // Density of drifting snow on virtual luv slope
			const double dL = -massErode / (EMS[nE-1].Rho);
			NDS[nE].z += dL;
			EMS[nE-1].L0 = EMS[nE-1].L = EMS[nE-1].L + dL;
			Xdata.cH += dL;
			NDS[nE].z += NDS[nE].u;
			NDS[nE].u = 0.0;
			NDS[nE].hoar = 0.;
			EMS[nE-1].M -= massErode;
			assert(EMS[nE-1].M>=0.); //mass must be positive
			Xdata.ErosionMass = massErode;
			forced_massErode = 0.;
		} else {
			Xdata.ErosionMass = 0.;
		}
		if (nErode > 0)
			Xdata.resize(nE);

		if (!alpine3d && SnowDrift::msg_erosion) { //messages on demand but not in Alpine3D
			if (Xdata.ErosionMass > 0.) {
				if (windward)
					prn_msg(__FILE__, __LINE__, "msg+", Mdata.date, "Eroding %d layer(s) w/ total mass %.3lf kg/m2 (windward: azi=%.1lf, slope=%.1lf)", nErode, Xdata.ErosionMass, Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());
				else
					prn_msg(__FILE__, __LINE__, "msg+", Mdata.date, "Eroding %d layer(s) w/ total mass %.3lf kg/m2 (azi=%.1lf, slope=%.1lf)", nErode, Xdata.ErosionMass, Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());
			}
		}
	// ... or, in case of no real erosion, check whether you can potentially erode at windward station using vw_drift.
	// This will as of 2023 contribute to the drift index VI24, to accommodate the ALPSolut use case!
	} else if ((snow_erosion == "TRUE") && (Xdata.ErosionLevel > Xdata.SoilNode) && Xdata.windward) {
		Sdata.drift = compMassFlux(EMS[Xdata.ErosionLevel], ustar_max, Xdata.meta.getSlopeAngle());  // kg m-1 s-1, main station, local vw && erosion level
		double virtuallyErodedMass = Sdata.drift * sn_dt / Hazard::typical_slope_length; // Convert to eroded snow mass in kg m-2
		if (virtuallyErodedMass > Constants::eps) {
			// Add (negative) value stored in Xdata.ErosionMass
			if (Xdata.ErosionMass < -Constants::eps)
				virtuallyErodedMass -= Xdata.ErosionMass;
			Sdata.mass[SurfaceFluxes::MS_WIND] = std::min(virtuallyErodedMass, EMS[Xdata.ErosionLevel].M); // use MS_WIND to carry virtually eroded mass
			// Now keep track of mass that either did or did not lead to erosion of full layer
			if (virtuallyErodedMass > EMS[Xdata.ErosionLevel].M) {
				virtuallyErodedMass -= EMS[Xdata.ErosionLevel].M;
				Xdata.ErosionLevel--;
			}
			Xdata.ErosionMass = -virtuallyErodedMass;
			Xdata.ErosionLevel = std::max(Xdata.SoilNode, std::min(Xdata.ErosionLevel, nE-1));
		} else {
			Xdata.ErosionMass = 0.;
		}
		if (!alpine3d && SnowDrift::msg_erosion) { //messages on demand but not in Alpine3D
			if (Xdata.ErosionLevel > nE-1)
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date, "Virtual erosion: ErosionLevel=%d did get messed up (nE-1=%d)", Xdata.ErosionLevel, nE-1);
		}
	} else {
		Xdata.ErosionMass = 0.;
	}
}


/**
 * @brief Erodes Elements from the top and computes the associated mass flux
 * Even so the code is quite obscure, it should cover all of the following cases:
 * -# Externally provided eroded mass (for example, by Alpine3D); see parameter forced_massErode
 * 
 * - SNOW_REDISTRIBUTION (SIMPLE/ADVANCED) erode at windward virtual slope (luv) and redistribute to lee slope (lee). 
 * 		Note that lee deposition is handled in Main.cc (setDataForCurrentTimeStep()). Will use VW_drift if available.
 * -# SNOW_EROSION (HS_DRIVEN, previously TRUE): using vw to erode the snow surface at the main station (flat field or slope).
 * -# SNOW_EROSION (REDEPOSIT/FREE): Erode on all aspects (except leeward if snow_redistribution is set to SIMPLE/ADVANCED) using vw_drift (if available and larger than vw) 
 * 		In case of REDEPOSIT the ErosionMass will be redeposited on the same slope (in Snowpack::Redepositsnow()).
 * -# SNOW_EROSION (VIRTUAL): using vw_drift (if available and larger than vw) to do virtual erosion either on flat field (#slopes == 1) or on windward slope (#slopes > 1).
 *
 * 	@note In case of SNOW_EROSION=HS_DRIVEN, erosion will also be controlled by mH and thus a measured snow depth (HS1) is required
 * 	@note If either measured snow depth is missing or the conditions for a real erosion are not fulfilled,
 *          the possibility of a virtual erosion will be considered using the ErosionLevel marker (virtual erodible layer).
 * @param Mdata
 * @param Xdata
 * @param Sdata
 * @param forced_massErode if greater than 0, force the eroded mass to the given value (instead of computing it)
*/
void SnowDrift::compSnowDrift(const CurrentMeteo& Mdata, SnowStation& Xdata, SurfaceFluxes& Sdata, double& forced_massErode) const
{
	// Supporting the deprecated old snow drift routine, for as long as necessary to verify correct working of new routine.
	if (snow_erosion == "TRUE" || snow_erosion == "FALSE") {
		compSnowDrift_deprecated(Mdata, Xdata, Sdata, forced_massErode);
		return;
	}

	size_t nE = Xdata.getNumberOfElements();
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;
	const bool no_snow = ((nE < Xdata.SoilNode+1) || (EMS[nE-1].theta[SOIL] > 0.));
	const bool no_wind_data = (Mdata.vw_drift == mio::IOUtils::nodata);
	Xdata.ErosionMass = 0.;

	if (no_snow || no_wind_data) {
		Xdata.ErosionMass = 0.;
		if (no_snow) {
			Xdata.ErosionLevel = Xdata.SoilNode;
			Sdata.drift = 0.;
		} else
			Sdata.drift = Constants::undefined;
		return;
	}

	// *******************************************************************
	// Check if conditions are met for erosion, based on snow_erosion and snow_redistribution, by station type (main or virtual slope).
	// If so, claculate mass flux of drifting snow.
	double massErode=0.; // Mass loss due to erosion
	const double ustar_max = (Mdata.vw>0.1) ? Mdata.ustar * Mdata.vw_drift / Mdata.vw : 0.; // Scale Mdata.ustar

	if (Xdata.sector == 0 ){ // Erode at main station if:
		if (fabs(forced_massErode) > Constants::eps2) {
			massErode = std::max(0., -forced_massErode); //negative mass is erosion
		// SLF operational setup: measured snow depth at main station, erode if discrepancy with calculated snow depth:
		} else if (enforce_measured_snow_heights && snow_erosion == "HS_DRIVEN" && (Xdata.mH > (Xdata.Ground + Constants::eps)) && ((Xdata.mH + 0.02) < Xdata.cH)) {
			Sdata.drift = compMassFlux(EMS[nE-1], Mdata.ustar, Xdata.meta.getSlopeAngle()); // kg m-1 s-1, main station, local vw
		} else if (snow_erosion == "FREE" || snow_erosion == "REDEPOSIT") {
			Sdata.drift = compMassFlux(EMS[nE-1], ustar_max, Xdata.meta.getSlopeAngle()); // kg m-1 s-1, main station, vw_drift 
		}
	} else if ( Xdata.sector > 0 ) { // Erode at virtual slope if:
		if (fabs(forced_massErode) > Constants::eps2) {
			massErode = std::max(0., -forced_massErode); //negative mass is erosion
		// in case of snow redistribution && erosion, erode everywhere except lee side. (luv-eroded snow will be deposited on the lee side in Main.cc).
		}else if ((snow_redistribution!="NONE") && (snow_erosion == "FREE" || snow_erosion == "REDEPOSIT")) {
			if (!Xdata.leeward) { 
				Sdata.drift = compMassFlux(EMS[nE-1], ustar_max, Xdata.meta.getSlopeAngle()); // kg m-1 s-1, vslope && vw_drift
			}
		} else if (snow_erosion == "FREE" || snow_erosion == "REDEPOSIT") { // when snow_redistribution is NONE, erode at ALL vslopes:
			Sdata.drift = compMassFlux(EMS[nE-1], ustar_max, Xdata.meta.getSlopeAngle()); // kg m-1 s-1, vslope && vw_drift  
		} else if ((snow_redistribution != "NONE") && (snow_erosion == "NONE" )){ //if snow_redistribution is true, but snow_erosion is not, erode only at luv (windward) side
			if (Xdata.windward ){
				Sdata.drift = compMassFlux(EMS[nE-1], ustar_max, Xdata.meta.getSlopeAngle()); // kg m-1 s-1, vslope && vw_drift
			}
		} // this should cover all cases of snow_erosion and snow_redistribution on virtual slopes. (except for SNOW_EROSION=VIRTUAL, see bottom of this function)
	}
	if ((Sdata.drift != Constants::undefined) && (Sdata.drift >0 )){ 
		//convert the drifing snow mass flux (kg m-1 s-1) to eroded snow mass ( in kg m-2 )
		massErode = Sdata.drift * sn_dt / Hazard::typical_slope_length; // Convert to eroded snow mass in kg m-2
	
		// *******************************************************************
		// Now erode the calculated mass from the top of the snowpack
		unsigned int nErode=0; // number of eroded elements
		// Check for limits to halt erosion. (t1au_thresh should limit erosion, but does a poor job. This allows for a hard rho value to halt erosion.) 
		if (erosion_limit != Constants::undefined && erosion_limit != 0.) {
			if (EMS[nE-1].Rho > erosion_limit) return;
		}
		if (massErode >= 0.95 * EMS[nE-1].M) {
			// Erode at most one element with a maximal error of +- 5 % on mass ...
			if ((snow_redistribution=="SIMPLE") && Xdata.windward) // in case of SIMPLE redistribution, use original density for deposition on virtual lee slope
				Xdata.rho_hn = EMS[nE-1].Rho;
			nE--;
			Xdata.cH -= EMS[nE].L;
			NDS[nE].hoar = 0.;
			Xdata.ErosionMass = EMS[nE].M;
			Xdata.ErosionLevel = std::min(nE-1, Xdata.ErosionLevel);
			nErode++;
			massErode -= EMS[nE].M;
			forced_massErode = -massErode;
		} else if (massErode > Constants::eps) { // ... or take away massErode from top element - partial real erosion
			if (fabs(EMS[nE-1].L * EMS[nE-1].Rho - EMS[nE-1].M) > 0.001) {
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date, "Inconsistent Mass:%lf   L*Rho:%lf", EMS[nE-1].M,EMS[nE-1].L*EMS[nE-1].Rho);
				EMS[nE-1].M = EMS[nE-1].L * EMS[nE-1].Rho;
				assert(EMS[nE-1].M>=0.); //mass must be positive
			}
			if ((snow_redistribution=="SIMPLE") && Xdata.windward) // in case of SIMPLE redistribution, use original density for deposition on virtual lee slope
				Xdata.rho_hn = EMS[nE-1].Rho; // Density of drifting snow on virtual luv slope
			const double dL = -massErode / (EMS[nE-1].Rho);
			NDS[nE].z += dL;
			EMS[nE-1].L0 = EMS[nE-1].L = EMS[nE-1].L + dL;
			Xdata.cH += dL;
			NDS[nE].z += NDS[nE].u;
			NDS[nE].u = 0.0;
			NDS[nE].hoar = 0.;
			EMS[nE-1].M -= massErode;
			assert(EMS[nE-1].M>=0.); //mass must be positive
			Xdata.ErosionMass = massErode;
			forced_massErode = 0.;
		} else {
			Xdata.ErosionMass = 0.;
		}
		if (nErode > 0)
			Xdata.resize(nE);

		if (!alpine3d && SnowDrift::msg_erosion) { //messages on demand but not in Alpine3D
			if (Xdata.ErosionMass > 0.) {
				if (Xdata.windward)
					prn_msg(__FILE__, __LINE__, "msg+", Mdata.date, "Eroding %d layer(s) w/ total mass %.3lf kg/m2 (windward: azi=%.1lf, slope=%.1lf)", nErode, Xdata.ErosionMass, Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());
				else
					prn_msg(__FILE__, __LINE__, "msg+", Mdata.date, "Eroding %d layer(s) w/ total mass %.3lf kg/m2 (azi=%.1lf, slope=%.1lf)", nErode, Xdata.ErosionMass, Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());
			}
		}
	}
	// ... or, in case of no real erosion, check whether you can potentially erode at windward station using vw_drift.
	// This will as of 2023 contribute to the drift index VI24, to accommodate the ALPSolut use case!
	if ((snow_erosion=="VIRTUAL") && (Xdata.ErosionLevel > Xdata.SoilNode) && Xdata.windward) {
		Sdata.drift = compMassFlux(EMS[Xdata.ErosionLevel], ustar_max, Xdata.meta.getSlopeAngle());  // kg m-1 s-1, main station, local vw && erosion level
		double virtuallyErodedMass = Sdata.drift * sn_dt / Hazard::typical_slope_length; // Convert to eroded snow mass in kg m-2
		if (virtuallyErodedMass > Constants::eps) {
			// Add (negative) value stored in Xdata.ErosionMass
			if (Xdata.ErosionMass < -Constants::eps)
				virtuallyErodedMass -= Xdata.ErosionMass;
			Sdata.mass[SurfaceFluxes::MS_WIND] = std::min(virtuallyErodedMass, EMS[Xdata.ErosionLevel].M); // use MS_WIND to carry virtually eroded mass
			// Now keep track of mass that either did or did not lead to erosion of full layer
			if (virtuallyErodedMass > EMS[Xdata.ErosionLevel].M) {
				virtuallyErodedMass -= EMS[Xdata.ErosionLevel].M;
				Xdata.ErosionLevel--;
			}
			Xdata.ErosionMass = -virtuallyErodedMass;
			Xdata.ErosionLevel = std::max(Xdata.SoilNode, std::min(Xdata.ErosionLevel, nE-1));
		} else {
			Xdata.ErosionMass = 0.;
		}
		if (!alpine3d && SnowDrift::msg_erosion) { //messages on demand but not in Alpine3D
			if (Xdata.ErosionLevel > nE-1)
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date, "Virtual erosion: ErosionLevel=%d did get messed up (nE-1=%d)", Xdata.ErosionLevel, nE-1);
		}
	}
}
