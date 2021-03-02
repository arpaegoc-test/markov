#include "../Distributions.h"
#include "FireGrowth.h"
#include "../TestResults.h"



FireGrowth::FireGrowth()
{

}

void FireGrowth::TestDistributions()
{
	Distribution *dist = 
				DistributionFactory::GetDistribution(DistributionFactory::Normal, 2.0, 5.0/2.0);

	double t= 10;
	double normalCdf = dist->cdf(t);

	double normalPdf = dist->pdf(t);

	// t = 10, mean = 2, sd=5
	std::cout << "Normal cdf = " << normalCdf << std::endl;	// 0.9452
	std::cout << "Normal pdf = " << normalPdf << std::endl;	// 0.0222



	Distribution *logdist = 
		DistributionFactory::GetDistribution(DistributionFactory::LogNormal, 5.18, 4.18/5.18);

	double lognormalCdf = logdist->cdf(t);
	double lognormalPdf = logdist->pdf(t);

	// t = 10, mean = 8.45, sd=0.78
	std::cout << "LogNormal cdf = " << lognormalCdf << std::endl;	// 9.695670493111511e-001
	std::cout << "LogNormal pdf = " << lognormalPdf << std::endl;	// 7.475127894932546e-002

	double mean = 5.18;
	double cov = 4.18/5.18;

	double lg = log(1+cov*cov);
	std::cout << log(mean)-lg/2 << std::endl;		//	2.12992
	std::cout << sqrt(lg) << std::endl;				//	0.092112

//1.39406
//0.708155


	return;

}

const double m21 = 2.0;		//Sustained to non-fire	-	exponential dist

const double m23 = 8.45;	//Sustained to vigorous	-	log-normal dist
const double c23 = 0.78 / m23;
const double m32 = 1;		//Vigorous to sustained	-	exponential

const double m34 = 5.55;	//Vigorous to interactive -	normal
const double c34 = 3.22 / m34;
const double m43 = 1.5;		//interactive to vigorous - exponential

const double m45 = 0.5;		//interactive to remote	-	exponential
const double m54 = 0.6;		//remote to interactive	-	exponential

const double m56 = 5.18;	//remote to fullroom	-	lognormal
const double c56 = 4.18 / m56;
const double m61 = 5.1;		//full-room to non-fire	-	exponential


void FireGrowth::RunModel()
{
	TestDistributions();
	//RunSMP();
}


void FireGrowth::RunSMP()
{
		std::ostringstream fpath;
		
		fpath << "FireGrowth\\FireGrowth_Markov.dat";

		int nStates = 6;

		JaggedMatrix *jmWbl = new JaggedMatrix(nStates);

		//Distribution *d21 = 
		//	DistributionFactory::GetDistribution(DistributionFactory::Weibull, m21,1);
		//jmWbl->AddDistribution(1,0,d21);


		//Distribution *d23 = 
		//	DistributionFactory::GetDistribution(DistributionFactory::LogNormal, m23,c23);
		//jmWbl->AddDistribution(1,2,d23);

		//
		//Distribution *d32 = 
		//	DistributionFactory::GetDistribution(DistributionFactory::Weibull, m32,1);
		//jmWbl->AddDistribution(2, 1, d32);


		//Distribution *d34 = 
		//	DistributionFactory::GetDistribution(DistributionFactory::Normal, m34,c34);
		//jmWbl->AddDistribution(2,3,d34);


		//Distribution *d43 = 
		//	DistributionFactory::GetDistribution(DistributionFactory::Weibull, m43,1);
		//jmWbl->AddDistribution(3,2,d43);


		//Distribution *d45 = 
		//	DistributionFactory::GetDistribution(DistributionFactory::Weibull, m45,1);
		//jmWbl->AddDistribution(3,4,d45);


		//Distribution *d54 = 
		//	DistributionFactory::GetDistribution(DistributionFactory::Weibull, m54,1);
		//jmWbl->AddDistribution(4,3,d54);


		//Distribution *d56 = 
		//	DistributionFactory::GetDistribution(DistributionFactory::LogNormal, m56,c56);
		//jmWbl->AddDistribution(4,5,d56);


		//Distribution *d61 = 
		//	DistributionFactory::GetDistribution(DistributionFactory::Weibull, m61,1);
		//jmWbl->AddDistribution(5,0,d61);



		//---EXPONENTIAL TEST-------------------------------------------------------------

		Distribution *d21 = 
			DistributionFactory::GetDistribution(DistributionFactory::Weibull, m21,1);
		jmWbl->AddDistribution(1,0,d21);


		Distribution *d23 = 
			DistributionFactory::GetDistribution(DistributionFactory::Weibull, m23,1);
		jmWbl->AddDistribution(1,2,d23);

		
		Distribution *d32 = 
			DistributionFactory::GetDistribution(DistributionFactory::Weibull, m32,1);
		jmWbl->AddDistribution(2, 1, d32);


		Distribution *d34 = 
			DistributionFactory::GetDistribution(DistributionFactory::Weibull, m34,1);
		jmWbl->AddDistribution(2,3,d34);


		Distribution *d43 = 
			DistributionFactory::GetDistribution(DistributionFactory::Weibull, m43,1);
		jmWbl->AddDistribution(3,2,d43);


		Distribution *d45 = 
			DistributionFactory::GetDistribution(DistributionFactory::Weibull, m45,1);
		jmWbl->AddDistribution(3,4,d45);


		Distribution *d54 = 
			DistributionFactory::GetDistribution(DistributionFactory::Weibull, m54,1);
		jmWbl->AddDistribution(4,3,d54);


		Distribution *d56 = 
			DistributionFactory::GetDistribution(DistributionFactory::Weibull, m56,1);
		jmWbl->AddDistribution(4,5,d56);


		Distribution *d61 = 
			DistributionFactory::GetDistribution(DistributionFactory::Weibull, m61,1);
		jmWbl->AddDistribution(5,0,d61);
//------------------------------------------------------------------------------------------



		jmWbl->Display();
		//jmWbl->DisplayInputs();

		SemiMarkovModel smpWbl(jmWbl);
		smpWbl.SetModelInput(45, 20000);
		smpWbl.SetupMatrices();

		smpWbl.ComputeStateProbabilities();

		vector<double> time = smpWbl.GetTimeVector();
		vector<double> p21 = smpWbl.GetStateProbability(1,0);
		vector<double> p22 = smpWbl.GetStateProbability(1,1);
		vector<double> p23 = smpWbl.GetStateProbability(1,2);
		vector<double> p24 = smpWbl.GetStateProbability(1,3);
		vector<double> p25 = smpWbl.GetStateProbability(1,4);
		vector<double> p26 = smpWbl.GetStateProbability(1,5);


		TestResults sysresults(fpath.str());
		
		sysresults.AddResult(time);
		sysresults.AddResult(p21);
		sysresults.AddResult(p22);
		sysresults.AddResult(p23);
		sysresults.AddResult(p24);
		sysresults.AddResult(p25);
		sysresults.AddResult(p26);

		sysresults.Serialize( );

}