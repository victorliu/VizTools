#define NOMINMAX
#define _USE_MATH_DEFINES
#include <float.h>
#include <limits>
#include <complex>
#include <cmath>
#include <cstdio>
#include <stdint.h>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "lodepng.h"
extern "C" {
#include "deltri.h"
#include "read_data_1d.h"
}
#include <tclap/CmdLine.h>

struct colormap_info_struct{
	double zmin, zmax;
	int type;
	// type = 0 means default grayscale
	// type = 1 means manually specified colors
	int ncolors;
	double *color;
} colormap_info;

void load_colormap(colormap_info_struct *info, const std::string &str){
	FILE *fp = fopen(str.c_str(), "rb");
	double *data = NULL;
	int ncols, nlines;
	int ret = read_data_1d(fp, &ncols, &nlines, &data);
	fclose(fp);

	info->ncolors = nlines;
	info->color = (double*)malloc(sizeof(double) * 5*info->ncolors);
	for(size_t i = 0; i < nlines; ++i){
		size_t j;
		for(j = 0; j < ncols; ++j){
			info->color[5*i+j] = data[i*ncols+j];
		}
		for( ; j < 5; ++j){
			info->color[5*i+j] = 1;
		}
	}
	double tfirst = info->color[5*(0)+0];
	double tlast = info->color[5*(info->ncolors-1)+0];
	double nrm = 1. / (tlast - tfirst);
	for(int i = 0; i < info->ncolors; ++i){
		info->color[5*i+0] = nrm * (info->color[5*i+0] - tfirst);
	}
	info->color[0] = 0;
	info->color[5*(info->ncolors-1)+0] = 1;

/*
	for(size_t i = 0; i < nlines; ++i){
		for(size_t j = 0; j < 5; ++j){
			fprintf(stderr, " %g", info->color[5*i+j]);
		}
		fprintf(stderr, "\n");
	}
*/
}
void load_colorscale(colormap_info_struct *info, const std::string &str){
	std::istringstream iss(str);
	std::vector<double> fields;
	std::copy(
		std::istream_iterator<double>(iss),
		std::istream_iterator<double>(),
		std::back_inserter<std::vector<double> >(fields)
	);
	info->ncolors = fields.size() / 5;
	info->color = (double*)malloc(sizeof(double) * 5*info->ncolors);
	for(size_t i = 0; i < fields.size(); i += 5){
		info->color[i+0] = fields[i+0];
		info->color[i+1] = fields[i+1];
		info->color[i+2] = fields[i+2];
		info->color[i+3] = fields[i+3];
		info->color[i+4] = fields[i+4];
	}
	double tfirst = info->color[5*(0)+0];
	double tlast = info->color[5*(info->ncolors-1)+0];
	double nrm = 1. / (tlast - tfirst);
	for(int i = 0; i < info->ncolors; ++i){
		info->color[5*i+0] = nrm * (info->color[5*i+0] - tfirst);
	}
	info->color[0] = 0;
	info->color[5*(info->ncolors-1)+0] = 1;
}

void colormap(const colormap_info_struct *info, double z, unsigned char *rgba){
	rgba[0] = 0;
	rgba[1] = 0;
	rgba[2] = 0;
	rgba[3] = 0;
	if(z == std::numeric_limits<double>::quiet_NaN()){ return; }
	if(!(z==z)){
		//fprintf(stderr, "nan in cm\n");
		return;
	}
	double t = (z - info->zmin) / (info->zmax - info->zmin);
	if(0 == info->type){
		rgba[0] = (int)(t*255+0.5);
		rgba[1] = rgba[0];
		rgba[2] = rgba[0];
		rgba[3] = 255;
	}else if(1 == info->type){
		for(int i = 0; i+1 < info->ncolors; ++i){
			//if(info->color[5*i+0] <= t && t <= info->color[5*(i+1)+0]){
			if(t <= info->color[5*(i+1)+0]){
				double tc = (t - info->color[5*i+0]) / (info->color[5*(i+1)+0] - info->color[5*i+0]);
				for(int p = 0; p < 4; ++p){
					double fc = (1.-tc) * info->color[5*i+1+p] + tc * info->color[5*(i+1)+1+p];
					rgba[p] = (int)(fc*255+0.5);
				}
				break;
			}
		}
	}
}


static int rand_int(int n) {
	const int limit = RAND_MAX - RAND_MAX % n;
	int rnd;
	do{ rnd = rand(); }while (rnd >= limit);
	return rnd % n;
}

/* go through and find some extremal points */
void fixup_data(int ndata, double *data, double scalex, double scaley){
	const double dirs[4][2] = { {1, 1}, {1, -1}, {-1, 1}, {-1, -1} };
	int ind[4] = { 0, 0, 0, 0 };
	double dot[4] = { -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX };

	/* Randomly permute */
	for(int i = ndata-1; i > 0; i--){
		const int k = rand_int(i+1);
		for(int j = 0; j < 3; ++j){
			double t = data[3*i+j];
			data[3*i+j] = data[3*k+j];
			data[3*k+j] = t;
		}
		// Find the points most in the direction of dirs
		for(int m = 0; m < 4; ++m){
			double cdot = dirs[m][0]*data[3*i+0] + dirs[m][1]*data[3*i+1];
			if(cdot > dot[m]){
				dot[m] = cdot;
				ind[m] = i;
			}
		}
	}
	std::sort(ind, ind+4);
	/*
	fprintf(stderr, "Extremal points:\n");
	for(int m = 0; m < 4; ++m){
		printf("%d\t%g, %g\n", ind[m], data[3*ind[m]+0], data[3*ind[m]+1]);
	}
	*/

	// swap them into place
	for(int m = 0; m < 4; ++m){
		if(m == ind[m]){ continue; }
		for(int j = 0; j < 3; ++j){
			double t = data[3*m+j];
			data[3*m+j] = data[3*ind[m]+j];
			data[3*ind[m]+j] = t;
		}
	}
	/*
	fprintf(stderr, "Extremal points:\n");
	for(int m = 0; m < 4; ++m){
		printf("%d\t%g, %g\n", m, data[3*m+0], data[3*m+1]);
	}
	*/

	for(int i = 0; i < ndata; ++i){
		data[3*i+0] *= scalex;
		data[3*i+1] *= scaley;
	}
}

double asinh(double x){
	// log(x+sqrt(1+x^2))
	if(x > 1){
		double y = 1./x;
		return log(x + fabs(x)*sqrt(1. + y*y));
	}else{
		return log(x + sqrt(x*x + 1.));
	}
}

class PositiveIntegerConstraint : public TCLAP::Constraint<int>{
	std::string description() const{
		return "Positive integer";
	}
	std::string shortID() const{
		return "pos-int";
	}
	bool check(const int &value) const{
		return value > 0;
	}
};

class NonNegativeDoubleConstraint : public TCLAP::Constraint<double>{
	std::string description() const{
		return "Non-negative real";
	}
	std::string shortID() const{
		return "non-neg-real";
	}
	bool check(const double &value) const{
		return value >= 0;
	}
};

int main(int argc, char *argv[]){
	double log_scale = 0;
	std::string in_filename;
	int xcol = 1, ycol = 2, zcol = 3;
	int oversample = 1;
	std::string colorscale;

	double xbounds[2] = { DBL_MAX, -DBL_MAX };
	double ybounds[2] = { DBL_MAX, -DBL_MAX };
	double zbounds[2] = { DBL_MAX, -DBL_MAX };

	double xrange, yrange, zrange;
	double azmax;

	int pxwidth, pxheight;

	double *data = NULL;
	int ndata;
	int cols[3];
	deltri DT = NULL;

	colormap_info_struct info;
	info.type = 0;
	info.color = NULL;

	try{
		PositiveIntegerConstraint positive_integer_constraint;
		NonNegativeDoubleConstraint nonneg_double_constraint;

		TCLAP::CmdLine cmd("Renders unorganized 2D point data", ' ', "1.0");

		TCLAP::ValueArg<int> arg_xcol("", "xcol", "Column containing x data", false, 1, &positive_integer_constraint);
		cmd.add(arg_xcol);
		TCLAP::ValueArg<int> arg_ycol("", "ycol", "Column containing y data", false, 2, &positive_integer_constraint);
		cmd.add(arg_ycol);
		TCLAP::ValueArg<int> arg_zcol("", "zcol", "Column containing z data", false, 3, &positive_integer_constraint);
		cmd.add(arg_zcol);

		TCLAP::ValueArg<double> arg_xmin("", "xmin", "Minimum value of x to use", false, 0, "float");
		cmd.add(arg_xmin);
		TCLAP::ValueArg<double> arg_xmax("", "xmax", "Minimum value of x to use", false, 1, "float");
		cmd.add(arg_xmax);
		TCLAP::ValueArg<double> arg_ymin("", "ymin", "Minimum value of y to use", false, 0, "float");
		cmd.add(arg_ymin);
		TCLAP::ValueArg<double> arg_ymax("", "ymax", "Minimum value of y to use", false, 1, "float");
		cmd.add(arg_ymax);
		TCLAP::ValueArg<double> arg_zmin("", "zmin", "Minimum value of z to use", false, 0, "float");
		cmd.add(arg_zmin);
		TCLAP::ValueArg<double> arg_zmax("", "zmax", "Minimum value of z to use", false, 1, "float");
		cmd.add(arg_zmax);

		TCLAP::ValueArg<int> arg_oversample("", "oversample", "Oversampling factor", false, 1, &positive_integer_constraint);
		cmd.add(arg_oversample);

		TCLAP::ValueArg<double> arg_pxmaxdim("", "maxdim-px", "Number of pixels to use for maximum dimension", false, 256, "integer");
		cmd.add(arg_pxmaxdim);
		TCLAP::ValueArg<double> arg_pxwidth("", "width-px", "Number of pixels to use for width (not including margins)", false, 0, "integer");
		cmd.add(arg_pxwidth);
		TCLAP::ValueArg<double> arg_pxheight("", "height-px", "Number of pixels to use for height (not including margins)", false, 0, "integer");
		cmd.add(arg_pxheight);

		TCLAP::ValueArg<std::string> arg_colorscale("", "colorscale", "Colorscale string", false, "", "sets of 5 numbers");
		cmd.add(arg_colorscale);
		TCLAP::ValueArg<std::string> arg_colormap("", "colormap", "Colormap file", false, "", "filename");
		cmd.add(arg_colormap);

		TCLAP::ValueArg<double> arg_log_scale("", "log-scale", "Output in logarithmic scale, value is the bias", false, 0, &nonneg_double_constraint);
		cmd.add(arg_log_scale);

		TCLAP::UnlabeledValueArg<std::string> arg_infile("input-file", "Input file", true, "", "Field output file");
		cmd.add(arg_infile);

		cmd.parse(argc, argv);

		in_filename = arg_infile.getValue();
		FILE *fp = fopen(in_filename.c_str(), "rb");
		cols[0] = arg_xcol.getValue() - 1;
		cols[1] = arg_ycol.getValue() - 1;
		cols[2] = arg_zcol.getValue() - 1;
		read_data_1d_cols(fp, 3, cols, &ndata, &data);
		fclose(fp);

		for(int i = 0; i < ndata; ++i){
			if(data[3*i+0] < xbounds[0]){ xbounds[0] = data[3*i+0]; }
			if(data[3*i+0] > xbounds[1]){ xbounds[1] = data[3*i+0]; }
			if(data[3*i+1] < ybounds[0]){ ybounds[0] = data[3*i+1]; }
			if(data[3*i+1] > ybounds[1]){ ybounds[1] = data[3*i+1]; }
			if(data[3*i+2] < zbounds[0]){ zbounds[0] = data[3*i+2]; }
			if(data[3*i+2] > zbounds[1]){ zbounds[1] = data[3*i+2]; }
		}
		xrange = xbounds[1] - xbounds[0];
		yrange = ybounds[1] - ybounds[0];
		zrange = zbounds[1] - zbounds[0];

		pxwidth = arg_pxwidth.getValue();
		pxheight = arg_pxheight.getValue();
		if(pxwidth > 0 && pxheight > 0){
			// use these specified dimensions
		}else{
			azmax = fabs(zbounds[1]);
			if(fabs(zbounds[0]) > azmax){ azmax = fabs(zbounds[0]); }
			if(xrange > yrange){
				pxwidth = arg_pxmaxdim.getValue();
				pxheight = (int)(yrange * (double)pxwidth/xrange + 0.5);
			}else{
				pxheight = arg_pxmaxdim.getValue();
				pxwidth = (int)(xrange * (double)pxheight/yrange + 0.5);
			}
		}

		log_scale = fabs(arg_log_scale.getValue());

		if(arg_colorscale.isSet()){
			info.type = 1;
			load_colorscale(&info, arg_colorscale.getValue());
		}
		if(arg_colormap.isSet()){
			info.type = 1;
			load_colormap(&info, arg_colormap.getValue());
		}

		info.zmin = zbounds[0];
		info.zmax = zbounds[1];
	}catch(TCLAP::ArgException &e){
		fprintf(stderr, "Error: %s for argument %s\n", e.error().c_str(), e.argId().c_str());
		return EXIT_FAILURE;
    }

	printf("Note: zrange = [ %.14g, %.14g ]\n", zbounds[0], zbounds[1]);

	const double dt_scale_x = (double)pxwidth/xrange;
	const double dt_scale_y = (double)pxheight/yrange;
	fixup_data(ndata, data, dt_scale_x, dt_scale_y);
	DT = deltri_new(ndata, data);

	// generate pixels
	unsigned char *img = (unsigned char*)malloc(sizeof(unsigned char) * 4*pxwidth*pxheight);
	const double dx = xrange / (double)(pxwidth * oversample);
	const double dy = yrange / (double)(pxheight * oversample);
	const double onrm = 1./ (double)(oversample*oversample);
	int hy = 0;
	for(int y = 0; y < pxheight; ++y){
		const double ty = ((double)(pxheight-1-y) + 0.5) / (double) pxheight;
		const double y0 = ybounds[0] + ty * yrange;
		for(int x = 0; x < pxwidth; ++x){
			const double tx = ((double)x + 0.5) / (double) pxwidth;
			const double x0 = xbounds[0] + tx * xrange;
			int h = hy;
			double value = 0;
			int nvals = 0;

			for(int oy = 0; oy < oversample; ++oy){
				const double toy = ((double)oy + 0.5) / (double) oversample - 0.5;
				for(int ox = 0; ox < oversample; ++ox){
					const double tox = ((double)ox + 0.5) / (double) oversample - 0.5;
					const double xy[2] = {
						dt_scale_x*(x0 + tox * dx),
						dt_scale_y*(y0 + toy * dy)
					};
					double intval;
					if(0 == DT_interpolate(DT, xy, &h, &intval)){
						nvals++;
						value += intval;
					}else{
					//	fprintf(stderr, "interp error, xy = %g, %g\n", xy[0], xy[1]);
					}
				}
			}
			if(0 == x){ hy = h; }
			if(nvals > 0){
				value /= (double)nvals;

				if(0 != log_scale){
					if(zbounds[0] > 0){
						value = std::log(value);
					}else{
						double bias = log_scale;
						double frac = std::abs(value)/azmax;
						double num = asinh(frac*bias);
						double denom = frac*asinh(bias);
						value *= num/denom;
					}
				}
			}else{
				value = std::numeric_limits<double>::quiet_NaN();
			}
			colormap(&info, value, &img[4*(x+y*pxwidth)+0]);
		}
	}
	DT_destroy(DT);

	{ // encode and save
		char *buffer = (char*)malloc(sizeof(char) * (32+strlen(in_filename.c_str())));
		sprintf(buffer, "%s.png", in_filename.c_str());

		LodePNG::Encoder encoder;
		std::vector<unsigned char> imgbuf;
		encoder.encode(imgbuf, img, pxwidth, pxheight);
		LodePNG::saveFile(imgbuf, buffer);

		free(buffer);
	}
	free(img);
	free(data);
	free(info.color);
    return EXIT_SUCCESS;
}
