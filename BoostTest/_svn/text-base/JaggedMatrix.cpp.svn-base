
#include "JaggedMatrix.h"
#include "RBDTools.h"

#include<cstdlib>

#include <boost/foreach.hpp>

//#include <boost/graph/graphviz.hpp>
//#include <boost/graph/adjacency_list.hpp>
//#inc/lude <boost/lexical_cast.hpp>
//#include <boost/graph/graph_traits.hpp> 

#include <algorithm> //transform
#include <numeric>
#include <functional>   // for multiplies

namespace ublas = boost::numeric::ublas;

//typedef boost::GraphvizDigraph Graph;
//typedef boost::graph_traits<boost::GraphvizDigraph>::edge_descriptor EdgeDescriptor;



using namespace boost::numeric::ublas;

JaggedMatrix::JaggedMatrix()
{
	_size = 0;
	_edgeCount = 0;
}


JaggedMatrix::JaggedMatrix(int n)
{
	_jaggedMatrix = matrix<boost::ptr_vector<Distribution>>(n,n);
	_meanMatrix = matrix<std::vector<double>>(n,n);
	_covMatrix = matrix<std::vector<double>>(n,n);
	_size = n;
	_edgeCount= 0;
    _adjacencyMatrix = matrix<int>(n,n);
}


JaggedMatrix::~JaggedMatrix()
{
	for(int i=0;i<_size;i++)
		for(int j=0;j<_size;j++)
		{
			if( _jaggedMatrix(i,j).size() == 0)
				continue;
			for(int c=0; c < _jaggedMatrix(i,j).size();c++)
				boost::ptr_vector<Distribution>::auto_type x= _jaggedMatrix(i,j).release(_jaggedMatrix(i,j).begin());
		}
}

void JaggedMatrix::AddDistribution(int r, int c, Distribution* distribution)
{
	if(_size == 0)
		throw std::exception("Null jagged matrix");

	_edgeCount++;
	_adjacencyMatrix(r,c) = 1;


	_jaggedMatrix(r,c).push_back(distribution);
}

void JaggedMatrix::AddDistribution(const TransitionInfo& info)
{
	AddDistribution(info._from, info._to, 
		DistributionFactory::GetDistribution(info._type, info._meantime, info._cov));
}

void JaggedMatrix::AddDistribution(int r, int c, 
								   double mean, double cov, 
								   DistributionFactory::DistributionType type)
{
	if(_size == 0)
		throw std::exception("Null jagged matrix");

	Distribution* distribution = 
		DistributionFactory::GetDistribution(type, mean, cov);

	AddDistribution(r,c,distribution);

	_meanMatrix(r,c).push_back(mean);
	_covMatrix(r,c).push_back(cov);


}



void JaggedMatrix::AddDistribution(int r, int c, 
								   double meantime, double cov, 
								   double repairtime, double repaircov,
								   DistributionFactory::DistributionType failuretype,
								   DistributionFactory::DistributionType repairtype
								   )
{
	if(_size == 0)
		throw std::exception("Null jagged matrix");

	AddDistribution(r,c, meantime, cov, failuretype);
	AddDistribution(r,c, repairtime, repaircov, repairtype);

}

void JaggedMatrix::DisplayInputs()
{
	for (int i = 0; i < _size; i++)
	{
		for (int j = 0; j < _size; j++)
		{
			if(_meanMatrix(i,j).size()==0)
				std::cout << "0\t";


			for(int c = 0; c < _meanMatrix(i,j).size(); c++)
				std::cout << 1/_meanMatrix(i,j)[c] << "\t";


		}
		std::cout << std::endl;
	}

	std::cout << std::endl << std::endl;


	for (int i = 0; i < _size; i++)
	{
		for (int j = 0; j < _size; j++)
		{
			if(_covMatrix(i,j).size()==0)
				std::cout << "1\t";

			for(int c = 0; c < _covMatrix(i,j).size(); c++)
				std::cout << _covMatrix(i,j)[c] << "\t";

			
		}
		std::cout << std::endl;
	}

}



// individual pdf function of each distribution
vector<double> JaggedMatrix::Cellpdf(double t, int r, int c)
{
	vector<double> result;

	if( _jaggedMatrix(r,c).empty())
		return result;

	int len = _jaggedMatrix(r, c).size();

	result = vector<double>(len);

	for (unsigned i = 0; i < _jaggedMatrix(r, c).size(); i++)
		result[i] = (_jaggedMatrix(r,c))[i].pdf(t);

	return result;
}


// individual pdf function of a row
vector<double> JaggedMatrix::pdf(double t, int r)
{
	vector<double> result(_size);
	std::fill(result.begin(), result.end(), 0.0);

	vector<double> cellpdf;

	for (unsigned i = 0; i < _size; i++)
	{
		cellpdf = Cellpdf(t, r, i);
		
		if(cellpdf.empty())
			continue;

		result[i] = cellpdf[0];
	}

	return result;
}


/// <summary>
/// individual Reliability function of each distribution
/// in cell [i,j]
/// </summary>
/// <param name="t">time step at which reliability is needed</param>
/// <param name="r">row of the jagged matrix</param>
/// <param name="c">row of the jagged matrix</param>
/// <returns>Reliability vector</returns>
vector<double> JaggedMatrix::CellReliability(double t, int r, int c)
{
	vector<double> result;

	if (_jaggedMatrix(r, c).empty())
		return result;

	int len = _jaggedMatrix(r, c).size();

	result = vector<double>(len);


	for( int i=0; i<len; i++)
		result[i] = 1 - _jaggedMatrix(r, c)[i].cdf(t);


	return result;
}

double JaggedMatrix::RowReliabilityProduct(double t, int r)
{
	vector<double> cellrel;
	double prod=1;
	for(int c=0; c<_size; c++)
	{
		cellrel = CellReliability(t,r,c);
		prod = std::accumulate(cellrel.begin(), cellrel.end(), prod, std::multiplies<double>());
	}
	return prod;
}

double JaggedMatrix::cdf(double t, int i, int j)
{
	if (_jaggedMatrix(i, j).empty())
		return 0;

	int len = _jaggedMatrix(i, j).size();

	if (len == 1)
		return _jaggedMatrix(i, j)[0].cdf(t);
}

/// <summary>
/// combined pdf function of a cell [r,c]
/// </summary>
/// <param name="t">time step at which pdf is needed</param>
/// <param name="r">row of the jagged matrix</param>
/// <param name="c">row of the jagged matrix</param>
/// <returns>pdf</returns>
double JaggedMatrix::pdf(double t, int i, int j)
{
	if (_jaggedMatrix(i, j).empty())
		return 0;

	int len = _jaggedMatrix(i, j).size();

	if (len == 1)
		return _jaggedMatrix(i, j)[0].pdf(t);

	vector<double> v1, v2;

	if (len == 2)
	{
		v1 = CellReliability(t, i, j);
		v2 = Cellpdf(t, i, j);

		//return v1[0]*v2[0] + v1[1]*v2[1];

		return v1[0]*v2[1] + v1[1]*v2[0];
	}

	if (len == 3)
	{
		v1 = CellReliability(t, i, j);
		v2 = Cellpdf(t, i, j);

		return v1[0] * v1[1] * v2[2] +
			v1[0] * v2[1] * v1[2] +
			v2[0] * v1[1] * v1[2];
	}

	throw std::exception("More than two compositions not supported for state reduction");
}


/// <summary>
/// pdf matrix keeping state reduction in mind.
/// a maximum of 3 distributions can be composed for
/// a single transition as of now
/// </summary>
/// <param name="t">timestep at which pdf is needed</param>
/// <returns>matrix of pdfs for each transition</returns>
matrix<double> JaggedMatrix::pdf(double t)
{
	matrix<double> pdfm = matrix<double>(_size, _size);
	double pdfval;

	for (int i = 0; i < _size; i++)
		for (int j = 0; j < _size; j++)
		{
			pdfval = pdf(t, i, j);
			pdfm(i, j) = pdfval;
		}
	return pdfm;
}


/// <summary>
/// Combined Reliability function of a single cell [i,j] array
/// </summary>
/// <param name="t">time step at which reliability is needed</param>
/// <param name="r">row of the jagged matrix</param>
/// <param name="c">row of the jagged matrix</param>
/// <returns>Reliability</returns>
double JaggedMatrix:: Reliability(double t, int r, int c)
{
	if (_jaggedMatrix(r, c).empty())
		return 1.0;

	int len = _jaggedMatrix(r, c).size();
	double result = 1.0;

	for (int i = 0; i < len; i++)
		result *= (1 - _jaggedMatrix(r, c)[i].cdf(t));

	return result;
}

/// <summary>
/// Combined Reliability function of a single row
/// </summary>
/// <param name="t">time step at which reliability is needed</param>
/// <param name="r">row of the jagged matrix</param>
/// <param name="c">row of the jagged matrix</param>
/// <returns>Reliability</returns>
vector<double> JaggedMatrix:: Reliability(double t, int r)
{
	vector<double> result(_size);

	for (int i = 0; i < _size; i++)
		result[i] = Reliability(t,r,i);

	return result;
}

/// <summary>
/// Reliability function of minimum of n random variables is
/// product of reliability functions of the n r.v.s
/// </summary>
/// <param name="t">time step at which reliability is needed</param>
/// <returns>Reliability matrix keeping state reduction in view</returns>
matrix<double> JaggedMatrix:: Reliability(double t)
{
	matrix<double> result = matrix<double>(_size, _size,1.0);
	double rel = 0.0;

	for (int i = 0; i < _size; i++)
		for (int j = 0; j < _size; j++)
		{
			if (_jaggedMatrix(i, j).empty())
			{
				result(i, j) = 1.0;
				continue;
			}

			rel = 1.0;
			BOOST_FOREACH(Distribution& dist, _jaggedMatrix(i,j))
				rel *= (1 - dist.cdf(t));
			result(i, j) = rel;
		}

		return result;
}

/// <summary>
/// CDF of minimum of n random variables is
/// CDF = 1 - RELIABILITY
/// </summary>
/// <param name="t">timestep at which cdf is needed</param>
/// <returns>cdf matrix ( = 1 - reliability)</returns>
matrix<double> JaggedMatrix::cdf(double t)
{
	matrix<double> rel = Reliability(t);
	matrix<double> result = matrix<double>(_size, _size);

	for(int i=0;i<_size;i++)
		for(int j=0;j<_size;j++)
			result(i,j) = 1 - rel(i,j);
	return result;
}


/// <summary>
/// Hazard or failure rate of all the distributions at time t
/// Simply finds pdf / rel
/// </summary>
/// <param name="t">Time t at which hazard is needed</param>
/// <returns>hazard at time t</returns>
matrix<double> JaggedMatrix:: Hazard(double t)
{
	matrix<double> pdfmatrix = pdf(t);
	matrix<double> rel = Reliability(t);

	matrix<double> hazard = element_div(pdfmatrix, rel);

	return hazard;
}

//same distribution for all the transitions
void JaggedMatrix::BuildModel(matrix<std::vector<double>> mean, 
										  matrix<std::vector<double>> cov,
										  DistributionFactory::DistributionType distributionType
										  )                                  
{
    int r, c, len, j;
    double vShape = 0.0;

	_jaggedMatrix = matrix<boost::ptr_vector<Distribution>>(_size, _size);
    _adjacencyMatrix = matrix<int>(_size, _size);


    for ( r=0; r< _size; r++)
        for(c=0; c<_size; c++)
        {
            if (mean(r, c).empty())
			{
				_adjacencyMatrix(r, c) = 0;
                continue;
			}

            len = mean(r, c).size();

            for (j = 0; j < len; j++)
            {
				AddDistribution(r,c, 
					 mean(r,c)[j], cov(r,c)[j], distributionType );
            }
        }
}

/// <summary>
/// Builds the SMP model with 
/// componentfailure and repair times, 
/// failure time covs are also supplied
/// </summary>
/// <param name="failureTime">mttf of each component</param>
/// <param name="covs">cov for each time to failure</param>
/// <param name="repairTime">repair time of each component</param>
JaggedMatrix::JaggedMatrix(
						   int componentCount,
						   double failureTime[], 
						   double failurecov[], 
						   double repairTime[],
						   double repaircov[],
						   DistributionFactory::DistributionType failuretype,
						   DistributionFactory::DistributionType repairtype
						   ) 
{

	std::vector<std::vector<double>> combos;

	RBDTools::Combinations(componentCount, combos);
	_size = combos.size();
	_jaggedMatrix = matrix<boost::ptr_vector<Distribution>>(_size, _size);

	_meanMatrix = matrix<std::vector<double>>(_size,_size);
	_covMatrix = matrix<std::vector<double>>(_size,_size);
    _adjacencyMatrix = matrix<int>(_size, _size);
	_edgeCount = 0;


    int i = 0, j = 0;
    int nonZeroIdx = 0;
	std::vector<double> target;
	double repcov=0;


    for (i = 0; i < _size; i++)
        for (j = 0; j < _size; j++)
        {
			if(RBDTools::IsMultipleTransitions(combos[i],combos[j]) || i==j)
	                continue;

            nonZeroIdx = 0;
            while (combos[i][nonZeroIdx] == combos[j][nonZeroIdx])
                nonZeroIdx++;


            //from 1 ->0 means a failure
            target = combos[j];
			
			if(j < i  )	//failure encountered
			{
				AddDistribution(i, j, 
								failureTime[nonZeroIdx], failurecov[nonZeroIdx], 
								   failuretype);
			}
			else 		//repair encountered
			{
				if(repairTime == NULL)
					continue;

				repcov = (repaircov == NULL) ? 1.0 : repaircov[nonZeroIdx];

				AddDistribution(i, j, 
					repairTime[nonZeroIdx], repcov, 
					repairtype);
			}

        }

}


JaggedMatrix::JaggedMatrix(matrix<std::vector<double>> mean,
						   matrix<std::vector<double>> cov,
						   DistributionFactory::DistributionType type) 
{
	_size = mean.size1();


	_meanMatrix = matrix<std::vector<double>>(_size,_size);
	_covMatrix = matrix<std::vector<double>>(_size,_size);


	BuildModel(mean, cov, type);
}

bool JaggedMatrix::PathExists(int fromState, int toState)
{
	if(_adjacencyMatrix(fromState, toState) == 1)
		return true;
	return false;
}

//time spent in inState on condition that next state is toState
double JaggedMatrix::SampleTimeSpent(int inState, int toState, double p)
{
	if(PathExists(inState, toState)==false)
		return 0;

	//inverting composed distributions is way too complex, just take 0th.
	return _jaggedMatrix(inState, toState)[0].InverseTransform();
}

void JaggedMatrix::Display()
{
	for(int i=0;i<_size;i++)
	{
		for(int j=0;j<_size;j++)
		{
			if( _jaggedMatrix(i,j).size() == 0)
				std::cout << "[0,0]\t";
			else
			{
				std::cout << "{";
				for(int c=0; c<_jaggedMatrix(i,j).size();c++)
					_jaggedMatrix(i,j)[c].DisplayParameters();
				std::cout << "}";
			}

		}

		std::cout << std::endl;
	}
}


/* upgraded boost from 1.40 to 1.44
*  This piece of code went broken after the upgrade
*  typedef boost::GraphvizDigraph doesn't work anymore
*/
void JaggedMatrix::Serialize(std::string filepath)
{
	/*
	std::vector<std::pair<int,int>> edges;
	std::vector<int> weights;
	

	using namespace boost;

	Graph g;

	graph_property<Graph, graph_graph_attribute_t>::
         type& graphAttr = get_property(g,graph_graph_attribute);

    graphAttr["rankdir"] = "RL";

	const boost::property_map<Graph, boost::vertex_attribute_t>::
         type& vertAttr = boost::get(boost::vertex_attribute, g);
    const boost::property_map<Graph, boost::edge_attribute_t>::
         type& edgeAttr = boost::get(boost::edge_attribute, g);


	vector<graph_traits<Graph>::vertex_descriptor> states(_size);

	

		for(int i=0; i < _size; i++)
		{
				states[i] = add_vertex(g);
				std::ostringstream out;
				out << i;
				vertAttr[states[i]]["label"] = out.str();
		}

    EdgeDescriptor edge;
	bool inserted;

	for(int i=0; i < _size; i++)
	{
		
			for(int j=0; j < _size; j++)
			{
				if( _jaggedMatrix(i,j).size()==0)
					continue;

				tie(edge, inserted) = add_edge(states[i], states[j], g);
				edgeAttr[edge]["label"] = _jaggedMatrix(i,j)[0].Notation();

			}
	}

 // Write the dot file
    boost::write_graphviz(filepath, g);
*/
}


matrix<double> JaggedMatrix::GetMarkovTransitionMatrix()
{
	matrix<double> result(_size, _size);

	for(int	r = 0; r < _size; r++)
	{
		for(int	c = 0; c < _size; c++)
		{
			result(r,c) = 0;

			if( _meanMatrix(r,c).size() != 0)
			{
				if( _meanMatrix(r,c)[0] != 0)
				{
					result(r, c) = 1 / _meanMatrix(r, c)[0];
				}
			}

			std::cout << result(r,c) << "\t";

		}
		std::cout << std::endl;

	}
	return result;
}