
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <fstream>
#include <sstream>

#include <stdlib.h>
#include <math.h>


/*
 *   Splits given string by given delimiter
 *   Returns substrings in vector
 *
 *   This function was inspired by code from stackoverflow.com (https://stackoverflow.com/a/14266139)
 *   It was posted under CC BY-SA 3.0 license (https://creativecommons.org/licenses/by-sa/3.0/)
 */

std::vector<std::string> split_by_delimiter(std::string str, char delim) 
{    
    std::vector<std::string> str_vector;
    size_t pos = 0;
    while ((pos = str.find(delim)) != std::string::npos) {
        str_vector.push_back(str.substr(0, pos));
        str.erase(0, pos + 1);
    }
    str_vector.push_back(str);
    
    return str_vector;
}

/*
 *   Returns squared double value
 */

inline static double sqr(double x) {
    return x * x;
}

/*
 *   Parses file and loads values to arrays
 */

void parse_csv(std::string file, std::vector<std::string>& quarters, std::vector<double>& inflation, std::vector<double>& unemployment) 
{
    std::ifstream data(file, std::ifstream::in);
    std::string line;
    std::vector<std::string> cells;
    // ignore header
    std::getline(data, line);
    // TODO: check if correct csv

    while(std::getline(data, line)) {
        cells = split_by_delimiter(line, ';');
        quarters.push_back(cells[0]);
        // TODO: handle exception
        inflation.push_back(stod(cells[1]));
        unemployment.push_back(stod(cells[2]));
    }
};

/*
 *   Finds and returns determinant of 3x3 matrix
 */

double determinant(double mat[3][3])
{
    double det;
    det = mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2])
          - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
          + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
    return det;
}

/*
 *   Finds the solution of system of 3 linear equations using cramer's rule
 */

void cramer(double quad[3][4], double& a, double& b, double& c)
{
    double d[3][3] = {
        { quad[0][0], quad[0][1], quad[0][2] },
        { quad[1][0], quad[1][1], quad[1][2] },
        { quad[2][0], quad[2][1], quad[2][2] },
    };
    
    double d1[3][3] = {
        { quad[0][3], quad[0][1], quad[0][2] },
        { quad[1][3], quad[1][1], quad[1][2] },
        { quad[2][3], quad[2][1], quad[2][2] },
    };
    
    double d2[3][3] = {
        { quad[0][0], quad[0][3], quad[0][2] },
        { quad[1][0], quad[1][3], quad[1][2] },
        { quad[2][0], quad[2][3], quad[2][2] },
    };
    
    double d3[3][3] = {
        { quad[0][0], quad[0][1], quad[0][3] },
        { quad[1][0], quad[1][1], quad[1][3] },
        { quad[2][0], quad[2][1], quad[2][3] },
    };
 
    double det = determinant(d);
    double det1 = determinant(d1);
    double det2 = determinant(d2);
    double det3 = determinant(d3);

    // apply Cramer's Rule
    a = det1 / det;
    b = det2 / det;
    c = det3 / det; 
}

/*
 *   Finds quadratic regression function of points in x and y arrays
 *   Sets values of coefficientss a, b, c and coefficients of determination r
 *   
 */

int quad_reg(const std::vector<double> x, const std::vector<double> y, double& c, double& b, double& a, double& r)
{
    int n = x.size();

    double sum_x = 0.0;
    double sum_x2 = 0.0;
    double sum_x3 = 0.0;
    double sum_x4 = 0.0;
    double sum_xy = 0.0;
    double sum_x2y = 0.0;
    double sum_y = 0.0;
    double sum_y2 = 0.0;
    double ssr = 0.0;
    double sst = 0.0;

    for (int i = 0; i < n; i++) { 
        sum_x  += x[i];       
        sum_x2 += sqr(x[i]);  
        sum_x3 += x[i] * sqr(x[i]);  
        sum_x4 += x[i] * x[i] * sqr(x[i]);  
        sum_xy += x[i] * y[i];
        sum_x2y += sqr(x[i]) * y[i];
        sum_y  += y[i];      
        sum_y2 += sqr(y[i]); 
    } 

    double quad[3][4] = {
        { (double)n, sum_x, sum_x2, sum_y },
        { sum_x,    sum_x2, sum_x3, sum_xy },
        { sum_x2,   sum_x3, sum_x4, sum_x2y },
    };
    
    cramer(quad, a, b, c);

    for (int i = 0; i < n; i++) { 
        ssr += sqr(y[i] - (a + b * x[i] + c * sqr(x[i])));
        sst += sqr(y[i] - (sum_y / n));
    } 
    r = 1 - (ssr / sst);

    return 0;
}

/*
 *   Finds linear regression function of points in x and y arrays
 *   Sets values of coefficients a, b and coefficients of determination r
 *   
 */

int lin_reg(const std::vector<double> x, const std::vector<double> y, double& a, double& b, double& r)
{
    int n = x.size();

    double sum_x = 0.0;
    double sum_x2 = 0.0;
    double sum_xy = 0.0;
    double sum_y = 0.0;
    double sum_y2 = 0.0;

    for (int i=0;i<n;i++){ 
        sum_x  += x[i];       
        sum_x2 += sqr(x[i]);  
        sum_xy += x[i] * y[i];
        sum_y  += y[i];      
        sum_y2 += sqr(y[i]); 
    } 

    double denom = (n * sum_x2 - sqr(sum_x));

    if (denom == 0) {
        // singular matrix
        a = 0;
        b = 0;
        r = 0;
        return 1;
    }

    a = (n * sum_xy  -  sum_x * sum_y) / denom;
    b = (sum_y * sum_x2  -  sum_x * sum_xy) / denom;
    r = (sum_xy - sum_x * sum_y / n) / sqrt((sum_x2 - sqr(sum_x)/n) * (sum_y2 - sqr(sum_y)/n));

    return 0; 
}

/*
 *   Finds phillips curves and generates script for gnuplot
 *   
 */

void generate_graphics(const std::string input, const std::vector<double> inflation, const std::vector<double> unemployment)
{
    double a, b, c, r;

    std::ofstream output("output.plt");

    output << "set terminal png size 500,500" << std::endl;
    output << "set key autotitle columnhead" << std::endl;
    output << "set datafile separator \";\"" << std::endl;
    output << "set xlabel \"Nezaměstanost [%]\"" << std::endl;
    output << "set ylabel \"Inflace [%]\"" << std::endl;
    

    lin_reg(unemployment, inflation, a, b, r);
    std::cout  << std::endl << "Linear regression:" << std::endl;
    std::cout << " y = (" << a << ")x + (" << b << ")" << std::endl;
    std::cout << " r^2: " << sqr(r) << std::endl << std::endl << std::endl;
    output << "f(x) = (" << a << ") * x + (" << b << ")" << std::endl;

    quad_reg(unemployment, inflation, a, b, c, r);
    std::cout << "Quadratic regression:" << std::endl;
    std::cout << " y = (" << a << ")x^2 + (" << b << ")x + (" << c  << ")" << std::endl;
    std::cout << " r^2: " << r << std::endl << std::endl;
    output << "g(x) = (" << a << ") * x * x + (" << b << ") * x + (" << c << ")" << std::endl;

    std::cout << "Script output.plt was created." << std::endl << "Run 'gnuplot output.plt' to generate graphics." << std::endl << std::endl;

    output << "set output 'linreg.png'" << std::endl;
    output << "set title 'Lineární regresní funkce'" << std::endl;
    output << "plot f(x) with lines, \"" << input << "\" using 3:2 notitle with points" << std::endl;


    output << "set output 'quadreg.png'" << std::endl;
    output << "set title 'Kvadratická regresní funkce'" << std::endl;
    output << "plot g(x) with lines, \"" << input << "\" using 3:2 notitle with points" << std::endl;

    output.close();                

    return;
}

/*
 *   Checks arguments
 *   
 */

int main(int argc, char *argv[])
{
    std::vector<std::string> quarters;
    std::vector<double> inflation, unemployment;

    if (argc != 2) {
        std::cout << "Error: invalid arguments" << std::endl << std::endl << "Program accepts exactly 1 argument (input csv file)." << std::endl;
        return EXIT_FAILURE;
    }

    parse_csv(argv[1], quarters, inflation, unemployment);

    generate_graphics(argv[1], inflation, unemployment);

    return EXIT_SUCCESS;
}

