#include "TestResults.h"

TestResults::TestResults(std::string filename)
	{
		_filename = filename;
	}


void TestResults::ClearResults()
	{
		_results.clear();
	}

void TestResults::AddResult(vector<double> column)
	{
		_results.push_back(column);
	}

void TestResults::Serialize()
	{
		std::ofstream ofs(_filename.c_str());

		int colCount = _results.size();
		int phiSize = _results[0].size();

		for( int i=0; i< phiSize; i++)
		{

			for(int j=0; j< colCount; j++)
				ofs << _results[j][i] << "\t";

			ofs << std::endl;
		}
		ofs.close();
	}