#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <vector>
#include <time.h>
#include <cmath>
#include <random>
#include <chrono>
#include "polar.h"


using namespace std;

int n = 10;
int info_length = 536;
int block_length = (1 << n);
int crc_length = 24;
int eff_info_length = info_length - crc_length;
//int frozen_length = block_length - eff_info_length;
int frozen_length = block_length - info_length;

int* bitreversal = new int[block_length];
int* information_location = new int[info_length];
int* frozen_location = new int[frozen_length];

int crc_poly[] = { 1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1 };
vector<int> crc(begin(crc_poly), end(crc_poly));
PolarCode polar(n, block_length, info_length, frozen_length, crc_length);

int sgn(double x)
{
	return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

double ptanh(double x)
{
	if (x <= -7.0)
		return -0.999998;
	else if (-7.0 < x && x <= -3.68)
		return 0.0004 * x - 0.9972;
	else if (-3.68 < x && x <= -1.82)
		return 0.0268 * x - 0.90;
	else if (-1.82 < x && x <= -1.24)
		return 0.1781 * x - 0.6247;
	else if (-1.24 < x && x <= -0.66)
		return 0.4605 * x - 0.2745;
	else if (-0.66 < x && x <= 0.66)
		return 0.8764 * x;
	else if (0.66 < x && x <= 1.24)
		return 0.4605 * x + 0.2745;
	else if (1.24 < x && x <= 1.82)
		return 0.1781 * x + 0.6247;
	else if (1.82 < x && x <= 3.68)
		return (0.0268 * x + 0.9);
	else if (3.68 < x && x <= 7.0)
		return 0.0004 * x + 0.9972;
	else
		return 0.999998;
}

double patanh(double x)
{
	if (-1 < x && x <= -0.999998)
		return -7.0;
	else if (-0.999998 < x && x <= -0.9987)
		return (x + 0.9972) / 0.0004;
	else if (-0.9987 < x && x <= -0.9488)
		return (x + 0.9) / 0.0268;
	else if (-0.9488 < x && x <= -0.8455)
		return (x + 0.6247) / 0.1781;
	else if (-0.8455 < x && x <= -0.5784)
		return (x + 0.2745) / 0.4605;
	else if (-0.5784 < x && x <= 0.5784)
		return x / 0.8764;
	else if (0.5784 < x && x <= 0.8455)
		return (x - 0.2745) / 0.4605;
	else if (0.8455 < x && x <= 0.9488)
		return (x - 0.6247) / 0.1781;
	else if (0.9488 < x && x <= 0.9987)
		return (x - 0.9) / 0.0268;
	else if (0.9987 < x && x <= 0.999998)
		return (x - 0.9972) / 0.0004;
	else
		return 7.0;

}

double gfunction(double L1, double L2)
{
	return 2.0 * patanh(ptanh(L1 / 2.0) * ptanh(L2 / 2.0));

	//    return  0.9375* sgn(L1) * sgn(L2) * min(fabs(L1),fabs(L2));

}

vector<int> bpsk(vector<int> code_bits)
{

	vector<int> bp(code_bits.size(), 0);

	for (int i = 0; i < code_bits.size(); i++)
	{
		bp.at(i) = 2 * code_bits.at(i) - 1;
	}

	return bp;
}

void readtxt(int block_length, int info_length, int frozen_length)
{
	fstream file;
	file.open("bitreversal.txt", ios::in);
	for (int i = 0; i < block_length; i++)
	{
		file >> bitreversal[i];
	}
	file.close();

	file.open("information_location.txt", ios::in);
	for (int i = 0; i < info_length; i++)
	{
		file >> information_location[i];
		//cout<< information_location[i]<<" ";
	}
	file.close();

	file.open("frozen_location.txt", ios::in);
	for (int i = 0; i < frozen_length; i++)
	{
		file >> frozen_location[i];
	}
	file.close();
}

vector<int> crcgen(vector<int> data)
{
	vector<int> adata(data.size() + crc.size() - 1, 0);

	for (int i = 0; i < data.size(); ++i)
	{
		adata.at(i) = data.at(i);
	}

	for (int i = 0; i < data.size(); ++i)
	{
		if (adata.at(i) == 1)
		{
			for (int j = i; j < i + crc.size(); ++j)
			{
				adata.at(j) = (adata.at(j) + crc.at(j - i)) % 2;
			}
		}
	}

	for (int i = 0; i < data.size(); ++i)
	{
		adata.at(i) = data.at(i);
	}

	return adata;
}

vector<int> PolarCode::encode(vector<int> block_bits)
{
	vector<int> info_bits_padded(_block_length, 0);
	vector <int> coded_bits(_block_length);

	for (int i = 0; i < _block_length; ++i)
	{
		info_bits_padded.at(i) = block_bits.at(i);
	}

	for (int iter = 0; iter < _n; ++iter)
	{
		int increment = (1 << iter);
		for (int j = 0; j < increment; j += 1)
		{
			for (int i = 0; i < _block_length; i += 2 * increment)
			{
				info_bits_padded.at(i + j) = ((info_bits_padded.at(i + j) + info_bits_padded.at(i + j + increment)) % 2);
			}
		}
	}


	for (int i = 0; i < _block_length; ++i)
	{
		coded_bits.at(i) = info_bits_padded.at(i);
	}

	return coded_bits;

}

bool PolarCode::crccheck(int* data)
{
	vector<int> crc_data(_info_length, 0);
	vector<int>  sys_decoded_bits(_info_length, 0);
	vector<int> sys_location_bits(_block_length, 0);
	vector<int> re_block_bits(_block_length, 0);
	vector<int> sys_code_bits(_block_length, 0);

	for (int i = 0; i < _info_length; ++i)
	{
		crc_data.at(i) = data[_channel_order_descending.at(i)];
	}

	for (int i = 0; i < _info_length; i++)
	{
		sys_decoded_bits.at(i) = crc_data.at(i);
	}

	for (int i = 0; i < _info_length; ++i)
	{
		sys_location_bits.at(information_location[i]) = sys_decoded_bits.at(i);
	}

	for (int i = 0; i < _block_length; ++i)
	{
		re_block_bits.at(i) = sys_location_bits.at(bitreversal[i]);
	}

	sys_code_bits = polar.encode(re_block_bits);
	//cout << "\nCRC chk\n ";
	for (int i = 0; i < info_length; i++)
	{
		crc_data.at(i) = sys_code_bits.at(bitreversal[information_location[i]]);
		//cout << crc_data.at(i) << " ";
	}

	int crc_size = _crc_length + 1;


	vector<int> check(crc_size - 1, 1);

	int pass = 0;
	bool crc_pass = true;

	for (int i = 0; i < _info_length - crc_size + 1; ++i)
	{
		if (crc_data[i] == 1)
		{
			for (int j = i; j < i + crc_size; ++j)
			{
				crc_data[j] = (crc_data[j] + crc.at(j - i)) % 2;
			}
		}
	}

	for (int i = 0; i < check.size(); ++i)
	{
		check.at(i) = crc_data[_info_length - crc_size + 1 + i];
	}

	for (int i = 0; i < check.size(); ++i)
	{
		pass = pass + check.at(i);
	}

	if (pass != 0)
	{
		crc_pass = false;
	}

	//cout << "\nPass : " << crc_pass << endl;
	return crc_pass;
}






vector<int> PolarCode::channel_order()
{
	_channel_order_descending.resize(_block_length);

	for (int i = 0; i < _block_length - _frozen_length; ++i)
	{
		_channel_order_descending.at(i) = information_location[i];
	}

	for (int j = _block_length - _frozen_length; j < _block_length; ++j)
	{
		_channel_order_descending.at(j) = frozen_location[j - (_block_length - _frozen_length)];
	}

	return _channel_order_descending;
}




void PolarCode::initializeDataStructures()
{
	while (_inactivePathIndices.size())   //new stack with capacity L
	{
		_inactivePathIndices.pop();
	};
	_activePath.resize(_list_size); //  new boolean array of size L

	_pathMetric_LLR.resize(_list_size);
	_arrayPointer_LLR.resize(_n + 1);
	for (int i = 0; i < _n + 1; ++i)
		_arrayPointer_LLR.at(i).resize(_list_size);

	_arrayPointer_C.resize(_n + 1);
	for (int i = 0; i < _n + 1; ++i) // new 2-D array of size (m+1) * L
		_arrayPointer_C.at(i).resize(_list_size);

	_arrayPointer_Info.resize(_list_size);

	_pathIndexToArrayIndex.resize(_n + 1);
	for (int i = 0; i < _n + 1; ++i)
		_pathIndexToArrayIndex.at(i).resize(_list_size);


	_inactiveArrayIndices.resize(_n + 1);
	for (int i = 0; i < _n + 1; ++i)
	{
		while (_inactiveArrayIndices.at(i).size())
		{
			_inactiveArrayIndices.at(i).pop();
		};
	}


	_arrayReferenceCount.resize(_n + 1);
	for (int i = 0; i < _n + 1; ++i)
		_arrayReferenceCount.at(i).resize(_list_size);

	for (int s = 0; s < _list_size; ++s)
	{
		_arrayPointer_Info.at(s) = new int[_block_length]();
		for (int lambda = 0; lambda < _n + 1; ++lambda)
		{
			_arrayPointer_LLR.at(lambda).at(s) = new double[(1 << (_n - lambda))]();
			_arrayPointer_C.at(lambda).at(s) = new int[2 * (1 << (_n - lambda))]();
			_arrayReferenceCount.at(lambda).at(s) = 0;
			_inactiveArrayIndices.at(lambda).push(s);
		}
	}


	for (int l = 0; l < _list_size; ++l)
	{
		_activePath.at(l) = 0;
		_inactivePathIndices.push(l);
		_pathMetric_LLR.at(l) = 0;
	}


}


int PolarCode::assignInitialPath()
{

	int  l = _inactivePathIndices.top();
	_inactivePathIndices.pop();
	_activePath.at(l) = 1;

	// Associate arrays with path index
	for (int lambda = 0; lambda < _n + 1; ++lambda)
	{
		int  s = _inactiveArrayIndices.at(lambda).top();
		_inactiveArrayIndices.at(lambda).pop();
		_pathIndexToArrayIndex.at(lambda).at(l) = s;
		_arrayReferenceCount.at(lambda).at(s) = 1;
	}
	return l;
}

int PolarCode::clonePath(int l)
{
	int l_p = _inactivePathIndices.top();
	_inactivePathIndices.pop();
	_activePath.at(l_p) = 1;

	_pathMetric_LLR.at(l_p) = _pathMetric_LLR.at(l);

	for (int lambda = 0; lambda < _n + 1; ++lambda)
	{
		int s = _pathIndexToArrayIndex.at(lambda).at(l);
		_pathIndexToArrayIndex.at(lambda).at(l_p) = s;
		_arrayReferenceCount.at(lambda).at(s)++;
	}
	return l_p;
}

void PolarCode::killPath(int l)
{
	_activePath.at(l) = 0;
	_inactivePathIndices.push(l);
	_pathMetric_LLR.at(l) = 0;

	for (int lambda = 0; lambda < _n + 1; ++lambda)
	{
		int s = _pathIndexToArrayIndex.at(lambda).at(l);
		_arrayReferenceCount.at(lambda).at(s)--;
		if (_arrayReferenceCount.at(lambda).at(s) == 0)
		{
			_inactiveArrayIndices.at(lambda).push(s);
		}
	}
}

double* PolarCode::getArrayPointer_LLR(int lambda, int  l)
{
	int  s = _pathIndexToArrayIndex.at(lambda).at(l);
	int s_p;
	if (_arrayReferenceCount.at(lambda).at(s) == 1)
	{
		s_p = s;
	}
	else
	{
		s_p = _inactiveArrayIndices.at(lambda).top();
		_inactiveArrayIndices.at(lambda).pop();

		//copy
		std::copy(_arrayPointer_C.at(lambda).at(s), _arrayPointer_C.at(lambda).at(s) + (1 << (_n - lambda + 1)), _arrayPointer_C.at(lambda).at(s_p));
		std::copy(_arrayPointer_LLR.at(lambda).at(s), _arrayPointer_LLR.at(lambda).at(s) + (1 << (_n - lambda)), _arrayPointer_LLR.at(lambda).at(s_p));

		_arrayReferenceCount.at(lambda).at(s)--;
		_arrayReferenceCount.at(lambda).at(s_p) = 1;
		_pathIndexToArrayIndex.at(lambda).at(l) = s_p;
	}
	return _arrayPointer_LLR.at(lambda).at(s_p);
}


int* PolarCode::getArrayPointer_C(int lambda, int  l)
{
	int  s = _pathIndexToArrayIndex.at(lambda).at(l);
	int s_p;
	if (_arrayReferenceCount.at(lambda).at(s) == 1)
	{
		s_p = s;
	}
	else
	{

		s_p = _inactiveArrayIndices.at(lambda).top();
		_inactiveArrayIndices.at(lambda).pop();

		//copy

		std::copy(_arrayPointer_LLR.at(lambda).at(s), _arrayPointer_LLR.at(lambda).at(s) + (1 << (_n - lambda)), _arrayPointer_LLR.at(lambda).at(s_p));
		std::copy(_arrayPointer_C.at(lambda).at(s), _arrayPointer_C.at(lambda).at(s) + (1 << (_n - lambda + 1)), _arrayPointer_C.at(lambda).at(s_p));

		_arrayReferenceCount.at(lambda).at(s)--;
		_arrayReferenceCount.at(lambda).at(s_p) = 1;
		_pathIndexToArrayIndex.at(lambda).at(l) = s_p;

	}
	return _arrayPointer_C.at(lambda).at(s_p);
}

void PolarCode::recursivelyCalcLLR(int lambda, int phi)
{
	if (lambda == 0)
		return;
	int psi = phi >> 1;
	if ((phi % 2) == 0)
		recursivelyCalcLLR(lambda - 1, psi);

	for (int l = 0; l < _list_size; ++l)
	{
		if (_activePath.at(l) == 0)
			continue;
		double* llr_lambda = getArrayPointer_LLR(lambda, l);
		double* llr_lambda_1 = getArrayPointer_LLR(lambda - 1, l);

		int* c_lambda = getArrayPointer_C(lambda, l);
		for (int beta = 0; beta < (1 << (_n - lambda)); ++beta)
		{
			if ((phi % 2) == 0)
			{
				//                if (40 > std::max(std::abs(llr_lambda_1[2 * beta]), std::abs(llr_lambda_1[2 * beta + 1])))
				//                {
				//                    llr_lambda[beta] = std::log ( (exp(llr_lambda_1[2 * beta] + llr_lambda_1[2 * beta + 1]) + 1) /
				//                                                  (exp(llr_lambda_1[2*beta]) + exp(llr_lambda_1[2*beta+1])));
				//                }
				//                else
				{
					llr_lambda[beta] = (double)((llr_lambda_1[2 * beta] < 0) ? -1 : (llr_lambda_1[2 * beta] > 0))*
						((llr_lambda_1[2 * beta + 1] < 0) ? -1 : (llr_lambda_1[2 * beta + 1] > 0))*
						std::min(std::fabs(llr_lambda_1[2 * beta]), std::fabs(llr_lambda_1[2 * beta + 1]));
				}
				//
				//                {
				//                llr_lambda[beta] = (double) gfunction(llr_lambda_1[2 * beta] ,llr_lambda_1[2 * beta + 1] );
				//                }
			}
			else
			{
				int  u_p = c_lambda[2 * beta];
				llr_lambda[beta] = (1 - 2 * u_p) * llr_lambda_1[2 * beta] + llr_lambda_1[2 * beta + 1];
			}

		}
	}
}


void PolarCode::recursivelyUpdateC(int lambda, int phi)
{

	int psi = phi >> 1;
	for (int l = 0; l < _list_size; ++l)
	{
		if (_activePath.at(l) == 0)
			continue;
		int* c_lambda = getArrayPointer_C(lambda, l);
		int* c_lambda_1 = getArrayPointer_C(lambda - 1, l);
		for (int beta = 0; beta < (1 << (_n - lambda)); ++beta)
		{
			c_lambda_1[2 * (2 * beta) + (psi % 2)] = (int)((c_lambda[2 * beta] + c_lambda[2 * beta + 1]) % 2);
			c_lambda_1[2 * (2 * beta + 1) + (psi % 2)] = c_lambda[2 * beta + 1];
		}
	}
	if ((psi % 2) == 1)
		recursivelyUpdateC((int)(lambda - 1), psi);

}

void PolarCode::continuePaths_FrozenBit(int phi)
{
	for (int l = 0; l < _list_size; ++l)
	{
		if (_activePath.at(l) == 0)
			continue;
		int* c_m = getArrayPointer_C(_n, l);

		c_m[(phi % 2)] = 0; // frozen value assumed to be zero

		double* llr_p = getArrayPointer_LLR(_n, l);

		_pathMetric_LLR.at(l) += log(1 + exp(-llr_p[0]));

		_arrayPointer_Info.at(l)[phi] = 0;
	}
}

void PolarCode::continuePaths_UnfrozenBit(int phi)
{

	std::vector<double>  probForks((unsigned long)(2 * _list_size));
	std::vector<double> probabilities;
	std::vector<int>  contForks((unsigned long)(2 * _list_size));


	int  i = 0;
	for (unsigned l = 0; l < _list_size; ++l)
	{
		if (_activePath.at(l) == 0)
		{
			probForks.at(2 * l) = NAN;
			probForks.at(2 * l + 1) = NAN;
		}
		else
		{

			double* llr_p = getArrayPointer_LLR(_n, l);
			probForks.at(2 * l) = -(_pathMetric_LLR.at(l) + log(1 + exp(-llr_p[0])));
			probForks.at(2 * l + 1) = -(_pathMetric_LLR.at(l) + log(1 + exp(llr_p[0])));

			probabilities.push_back(probForks.at(2 * l));
			probabilities.push_back(probForks.at(2 * l + 1));

			i++;
		}
	}

	int  rho = _list_size;
	if ((2 * i) < _list_size)
		rho = (int)2 * i;

	for (int l = 0; l < 2 * _list_size; ++l)
	{
		contForks.at(l) = 0;
	}
	std::sort(probabilities.begin(), probabilities.end(), std::greater<double>());

	double threshold = probabilities.at((unsigned long)(rho - 1));
	int num_paths_continued = 0;

	for (int l = 0; l < 2 * _list_size; ++l)
	{
		if (probForks.at(l) > threshold)
		{
			contForks.at(l) = 1;
			num_paths_continued++;
		}
		if (num_paths_continued == rho)
		{
			break;
		}
	}

	if (num_paths_continued < rho)
	{
		for (int l = 0; l < 2 * _list_size; ++l)
		{
			if (probForks.at(l) == threshold)
			{
				contForks.at(l) = 1;
				num_paths_continued++;
			}
			if (num_paths_continued == rho)
			{
				break;
			}
		}
	}

	for (unsigned l = 0; l < _list_size; ++l)
	{
		if (_activePath.at(l) == 0)
			continue;
		if (contForks.at(2 * l) == 0 && contForks.at(2 * l + 1) == 0)
			killPath(l);
	}

	for (unsigned l = 0; l < _list_size; ++l)
	{
		if (contForks.at(2 * l) == 0 && contForks.at(2 * l + 1) == 0)
			continue;
		int* c_m = getArrayPointer_C(_n, l);

		if (contForks.at(2 * l) == 1 && contForks.at(2 * l + 1) == 1)
		{

			c_m[(phi % 2)] = 0;
			int l_p = clonePath(l);
			c_m = getArrayPointer_C(_n, l_p);
			c_m[(phi % 2)] = 1;

			std::copy(_arrayPointer_Info.at(l), _arrayPointer_Info.at(l) + phi, _arrayPointer_Info.at(l_p));
			_arrayPointer_Info.at(l)[phi] = 0;
			_arrayPointer_Info.at(l_p)[phi] = 1;


			double* llr_p = getArrayPointer_LLR(_n, l);
			_pathMetric_LLR.at(l) += log(1 + exp(-llr_p[0]));
			llr_p = getArrayPointer_LLR(_n, l_p);
			_pathMetric_LLR.at(l_p) += log(1 + exp(llr_p[0]));


		}
		else
		{
			if (contForks.at(2 * l) == 1)
			{
				c_m[(phi % 2)] = 0;
				_arrayPointer_Info.at(l)[phi] = 0;


				double* llr_p = getArrayPointer_LLR(_n, l);
				_pathMetric_LLR.at(l) += log(1 + exp(-llr_p[0]));

			}
			else
			{
				c_m[(phi % 2)] = 1;
				_arrayPointer_Info.at(l)[phi] = 1;

				double* llr_p = getArrayPointer_LLR(_n, l);
				_pathMetric_LLR.at(l) += log(1 + exp(llr_p[0]));

			}
		}
	}

}


std::vector<int> PolarCode::decode_scl_llr(std::vector<double> llr, int list_size)
{

	_list_size = list_size;


	_channel_order_descending = channel_order();

	frozen();


	initializeDataStructures();

	int  l = assignInitialPath();

	double* llr_0 = getArrayPointer_LLR(0, l);

	for (int beta = 0; beta < _block_length; ++beta)
	{
		llr_0[beta] = llr.at(beta);
	}

	for (int phi = 0; phi < _block_length; ++phi)
	{


		recursivelyCalcLLR(_n, phi);



		if (_frozen_bits.at(phi) == 1)
			continuePaths_FrozenBit(phi);
		else
			continuePaths_UnfrozenBit(phi);

		if ((phi % 2) == 1)
			recursivelyUpdateC(_n, phi);

	}

	l = findMostProbablePath((bool)_crc_length);
	//cout << "\nAns path : " << l<<endl;
	int* c_0 = _arrayPointer_Info.at(l);
	std::vector<int> deocded_info_bits(_info_length);
	for (int beta = 0; beta < _info_length; ++beta)
		deocded_info_bits.at(beta) = c_0[_channel_order_descending.at(beta)];

	for (int s = 0; s < _list_size; ++s)
	{
		delete[] _arrayPointer_Info.at(s);
		for (int lambda = 0; lambda < _n + 1; ++lambda)
		{

			delete[] _arrayPointer_LLR.at(lambda).at(s);
			delete[] _arrayPointer_C.at(lambda).at(s);
		}
	}

	return deocded_info_bits;
}



int PolarCode::findMostProbablePath(bool check_crc)
{

	int  l_p = 0;
	double p_llr = numeric_limits<double>::max();
	bool path_with_crc_pass = false;

	for (int l = 0; l < _list_size; ++l)
	{

		if (_activePath.at(l) == 0)
			continue;

		if ((check_crc) && (!crccheck(_arrayPointer_Info.at(l))))
			continue;

		path_with_crc_pass = true;



		if (_pathMetric_LLR.at(l) < p_llr)
		{
			p_llr = _pathMetric_LLR.at(l);
			l_p = l;
		}


	}

	if (path_with_crc_pass)
		return l_p;
	else
		return findMostProbablePath(false);

}

void PolarCode::frozen()
{
	_frozen_bits.resize(_block_length);

	for (int i = 0; i < _frozen_length; i++)
	{
		_frozen_bits.at(frozen_location[i]) = 1;
	}

}






int main()
{

	srand((unsigned)time(NULL) + rand());
	std::mt19937 generator(static_cast<long unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
	std::normal_distribution<double> unif(0.0, 1.0);

	
	int list_size = 2;
	cout << "Systematic SCL + CRC\nList = ";
	cin >> list_size;
	double SNR;
	cout << "SNR = ";
	cin >> SNR;
	vector<int> info_bits;
	vector<int> re_info_bits(info_length, 0);
	vector<int> location_bits(block_length, 0);
	vector<int> re_location_bits(block_length, 0);
	vector<int> sys_location_bits(block_length, 0);
	vector<int> block_bits(block_length, 0);
	vector<int> re_block_bits(block_length, 0);
	vector<int> code_bits(block_length, 0);
	vector<int> sys_code_bits(block_length, 0);
	vector<int> bpsk_bits(block_length, 0);
	vector<double> receive_bits(block_length, 0);
	vector<double> receive_llr(block_length, 0);
	vector<int> decode_bits;
	vector<int>  sys_decoded_bits(info_length, 0);
	vector<int> crc_data(eff_info_length, 0);
	readtxt(block_length, info_length, frozen_length);

	double noise;
	double sigma;
	double coderate = (double)eff_info_length / block_length;
	double countimes = 35000000;
	int errorcount;
	double errorate;
	int blockerrorcount;
	double blockerrorate;

	//for (double SNR = 1; SNR <= 3; SNR += 0.5)
	{
		sigma = pow(1.0 / (2.0 * coderate * pow(10.0, SNR / 10.0)), 0.5);

		errorcount = 0;
		blockerrorcount = 0;
		double early_times;
		for (int times = 0; times < countimes; times++)
		{

			for (int i = 0; i < crc_data.size(); i++)
			{
				crc_data.at(i) = 2 * rand() / (RAND_MAX + 1.0);

			}

			info_bits = crcgen(crc_data);
			//for (int i = 0; i < info_bits.size(); ++i)
			//{
			//	cout << info_bits.at(i) << " ";
			//}
			for (int i = 0; i < info_length; ++i)
			{
				location_bits.at(information_location[i]) = info_bits.at(i);
			}

			for (int i = 0; i < block_length; ++i)
			{
				block_bits.at(i) = location_bits.at(bitreversal[i]);
			}

			code_bits = polar.encode(block_bits);

			for (int i = 0; i < info_length; i++)
			{
				re_info_bits.at(i) = code_bits.at(bitreversal[information_location[i]]);
			}

			for (int i = 0; i < info_length; ++i)
			{
				re_location_bits.at(information_location[i]) = re_info_bits.at(i);
			}

			for (int i = 0; i < block_length; ++i)
			{
				re_block_bits.at(i) = re_location_bits.at(bitreversal[i]);
			}

			code_bits = polar.encode(re_block_bits);

			bpsk_bits = bpsk(code_bits);



			for (int i = 0; i < block_length; ++i)
			{
				noise = sigma * unif(generator);
				receive_bits.at(i) = bpsk_bits.at(i) + noise;
			}

			for (int i = 0; i < receive_bits.size(); ++i)
			{
				receive_llr.at(i) = (-2.0 * receive_bits.at(i)) / (sigma * sigma);
			}

			decode_bits = polar.decode_scl_llr(receive_llr, list_size);
			//cout << "\ndecoded\n";
			for (int i = 0; i < info_length; i++)
			{
				
				sys_decoded_bits.at(i) = decode_bits.at(i);
			}

			for (int i = 0; i < info_length; ++i)
			{
				sys_location_bits.at(information_location[i]) = sys_decoded_bits.at(i);
			}

			for (int i = 0; i < block_length; ++i)
			{
				re_block_bits.at(i) = sys_location_bits.at(bitreversal[i]);
			}

			sys_code_bits = polar.encode(re_block_bits);

			for (int i = 0; i < info_length; i++)
			{
				decode_bits.at(i) = sys_code_bits.at(bitreversal[information_location[i]]);
				//cout << decode_bits.at(i) << " ";
			}


			for (int i = 0; i < info_length; i++)
			{
				if (decode_bits.at(i) != info_bits.at(i))
				{
					errorcount = errorcount + 1;
				}
			}

			for (int i = 0; i < info_length; i++)
			{
				if (decode_bits.at(i) != info_bits.at(i))
				{
					blockerrorcount = blockerrorcount + 1;
					break;
				}
			}

			early_times = (double)times;
			if (early_times >= 1) {
				cout << "\rSNR : " << SNR << " List : " << list_size << " err : " << errorcount << " ber : " << (double)(errorcount / (double)(times * eff_info_length * 1))
					<< " block err : " << blockerrorcount << " bler : " << (double)(blockerrorcount / ((double)times * 1.0))<<" Times : "<< early_times <<"    ";
			}
			if (blockerrorcount >= 500) {
				break;
			}

		}


		errorate = (double)(errorcount / (early_times * eff_info_length * 1.0));
		blockerrorate = (double)(blockerrorcount / (early_times * 1.0));

		cout << endl;
		cout << "SNR: " << SNR << endl;
		cout << "errorcount: " << errorcount << endl;
		cout << "error rate: " << errorate << endl;
		cout << "block errorcount: " << blockerrorcount << endl;
		cout << "block errorate: " << blockerrorate << endl;
		cout << endl;

	}


	system("pause");


	return 0;
}
