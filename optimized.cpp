#include <iostream>
#include <fstream>
#include <limits>
#include <iomanip>
#include <string>

using namespace std;

const double EPSILON = 1e-8;

const int DEC_PRECISION = 3;

int main(int argc, char* argv[]) {

    // Read input file
    ifstream file(argv[1]);

    // Read sizes
    int num_rest, num_vars;
    file >> num_rest >> num_vars;

    int total_cols = 2 * num_rest + num_vars + 1;

    // Create tableau
    double** tableau = new double*[num_rest + 2];
    for (int i = 0; i < num_rest + 2; ++i) {
        tableau[i] = new double[2 * num_rest + num_vars + 1];
    }

    // Filling the 1st row of the tableau:
    // Fill the 2 phase and certificate first row with 0's
    for (int i = 0; i < 2 * num_rest; ++i) {
        tableau[0][i] = 0;
    }
    // Read cost func and fill the tableau
    double temp;
    for (int i = 0; i < num_vars; ++i) {
        file >> temp;
        tableau[0][2 * num_rest + i] = -temp;
    }
    // Fill the top right element of the tableau with 0
    tableau[0][total_cols - 1] = 0;

    // Filling the middle of the tableau
    // Read restrictions
    for (int row = 0; row < num_rest; ++row) {
        for (int col = 0; col < num_vars + 1; ++col) {
            file >> temp;
            tableau[row + 1][2 * num_rest + col] = temp;
        }
    }
    // Fill the 2 phase and certificate identity matrix 
    for (int row = 0; row < num_rest; ++row) {
        for (int col = 0; col < num_rest; ++col) {
            if (row == col) {
                tableau[row + 1][col] = 1;
                if (tableau[row + 1][total_cols - 1] < 0) {
                    // case where b is negative
                    tableau[row + 1][col + num_rest] = -1;
                } else {
                    tableau[row + 1][col + num_rest] = 1;
                }
            } else {
                tableau[row + 1][col] = 0;
                tableau[row + 1][col + num_rest] = 0;
            }
        }
    }

    // Filling the last line of the tableau (for the 2 phase)
    // Fill the 2 phase and certificate with 1's
    for (int i = 0; i < num_rest; ++i) {
        tableau[num_rest + 1][i] = 1;
    }
    for (int i = 0; i < num_rest; ++i) {
        tableau[num_rest + 1][num_rest + i] = 0;
    }
    // Fill the cost func of the last row
    for (int i = 0; i < num_vars; ++i) {
        tableau[num_rest + 1][2 * num_rest + i] = 0;
    }
    // Fill the bottom right element of the tableau with 0
    tableau[num_rest + 1][total_cols - 1] = 0;
    file.close();

    // Make the negative b restrictions positive
    for (int row = 0; row < num_rest; ++row) {
        if (tableau[row + 1][total_cols - 1] < 0) {
            for (int col = 0; col < num_vars + 1; ++col) {
                if (tableau[row + 1][2 * num_rest + col] != 0) {
                    tableau[row + 1][2 * num_rest + col] = -tableau[row + 1][2 * num_rest + col];
                }
            }
        }
    }

    // Create basis
    int* basis = new int[num_rest];
    // Fill the basis with the columns of the 1st phase identity matrix
    for (int i = 0; i < num_rest; ++i) {
        basis[i] = i;
    }

    string input_filename(argv[1]);
    size_t start_pos = input_filename.find_last_of("/") + 1;
    size_t end_pos = input_filename.find_last_of(".");
    string base_name = input_filename.substr(start_pos, end_pos - start_pos);

    string output_filename = base_name + "_saida.txt";

    // Create the full path for the output file in the root directory
    string output_dir = "./" + output_filename;

    ofstream output(output_dir);
    output << fixed << setprecision(DEC_PRECISION);

    // 1st phase of the simplex
    // Eliminate 1's from the 1st phase
    for (int row = 0; row < num_rest; ++row) {
        for (int col = 0; col < total_cols; ++col) {
            tableau[num_rest + 1][col] -= tableau[row + 1][col];
        }
    }
    // Simplex 1st phase
    while (true) {
        // Find the pivot column
        int col_pivot = -1;
        double lowest_val = numeric_limits<double>::infinity();
        for (int i = 0; i < num_vars; ++i) {
            if (tableau[num_rest + 1][2 * num_rest + i] < -EPSILON && tableau[num_rest + 1][2 * num_rest + i] < (lowest_val - EPSILON)) {
                col_pivot = 2 * num_rest + i;
                lowest_val = tableau[num_rest + 1][2 * num_rest + i];
            }
        }
        if (col_pivot == -1) {
            break;
        }
        // Find the pivot row
        int row_pivot = -1;
        double min_ratio = numeric_limits<double>::infinity();
        double ratio = 0;
        for (int i = 0; i < num_rest; ++i) {
            if (tableau[i + 1][col_pivot] > EPSILON) {
                ratio = tableau[i + 1][total_cols - 1] / tableau[i + 1][col_pivot];
                if ((ratio < (min_ratio - EPSILON) || ((fabs(ratio - min_ratio) < EPSILON) && (basis[i] < basis[row_pivot - 1]))) && ratio >= -EPSILON) {
                    min_ratio = ratio;
                    row_pivot = i + 1;
                }
            }
        }
        if (row_pivot == -1) {
            break;
        }
        // Normalize the pivot row
        double norm_factor = tableau[row_pivot][col_pivot];
        for (int i = 0; i < total_cols; ++i) {
            tableau[row_pivot][i] /= norm_factor;
        }
        // Zero out the pivot column
        double factor = 0;
        for (int row = 0; row < num_rest + 2; ++row) {
            if (row != row_pivot) {
                factor = tableau[row][col_pivot];
                if (fabs(factor) > EPSILON) {
                    for (int col = 0; col < total_cols; ++col) {
                        tableau[row][col] -= factor * tableau[row_pivot][col];
                    }
                }
            }
        }
        // Update the basis
        basis[row_pivot - 1] = col_pivot;
    }
    // Check if the 1st phase is infeasible
    if (fabs(tableau[num_rest + 1][total_cols - 1]) > EPSILON) {
        //output << fixed << setprecision(DEC_PRECISION);
        output << "inviavel" << endl;
        // Write certificate of unfeasability
        for (int i = 0; i < num_rest; ++i) {
            output << tableau[num_rest + 1][num_rest + i] << " ";
        }
        output << endl;


    } else {
        // Simplex 2nd phase
        bool is_unbounded = false;
        int col_pivot = -1;
        while (true) {
            // Find the pivot column
            col_pivot = -1;
            double lowest_val = numeric_limits<double>::infinity();
            for (int i = 0; i < num_vars; ++i) {
                if (tableau[0][2 * num_rest + i] < -EPSILON && tableau[0][2 * num_rest + i] < (lowest_val - EPSILON)) {
                    col_pivot = 2 * num_rest + i;
                    lowest_val = tableau[0][2 * num_rest + i];
                }
            }
            if (col_pivot == -1) {
                break;
            }
            // Find the pivot row
            int row_pivot = -1;
            double min_ratio = numeric_limits<double>::infinity();
            double ratio = 0;
            for (int i = 0; i < num_rest; ++i) {
                if (tableau[i + 1][col_pivot] > EPSILON) {
                    ratio = tableau[i + 1][total_cols - 1] / tableau[i + 1][col_pivot];
                    if ((ratio < (min_ratio - EPSILON) || ((fabs(ratio - min_ratio) < EPSILON) && (basis[i] < basis[row_pivot - 1]))) && ratio >= -EPSILON) {
                        min_ratio = ratio;
                        row_pivot = i + 1;
                    }
                }
            }
            if (row_pivot == -1) {
                is_unbounded = true;
                break;
            }
            // Normalize the pivot row
            double norm_factor = tableau[row_pivot][col_pivot];
            for (int i = 0; i < total_cols; ++i) {
                tableau[row_pivot][i] /= norm_factor;
            }
            // Zero out the pivot column
            double factor = 0;
            for (int row = 0; row < num_rest + 2; ++row) {
                if (row != row_pivot) {
                    factor = tableau[row][col_pivot];
                    if (fabs(factor) > EPSILON) {
                        for (int col = 0; col < total_cols; ++col) {
                            tableau[row][col] -= factor * tableau[row_pivot][col];
                        }
                    }
                }
            }
            // Update the basis
            basis[row_pivot - 1] = col_pivot;
        }
        if (!is_unbounded) {
            output << "otimo" << endl;
            //output << fixed << setprecision(DEC_PRECISION);
            // Write the optimal value
            output << tableau[0][total_cols - 1] << endl;
            // Write the optimal solution
            int basic_row_index = -1;
            for (int i = 0; i < num_vars; ++i) {
                basic_row_index = -1;
                for (int j = 0; j < num_rest; ++j) {
                    if (basis[j] == 2 * num_rest + i) {
                        basic_row_index = j;
                        break;
                    }
                }
                if (basic_row_index != -1) {
                    output << tableau[basic_row_index + 1][total_cols - 1] << " ";
                } else {
                    output << "0.000 ";
                }
            }
            output << endl;
            // Write certificate of optimality
            for (int i = 0; i < num_rest; ++i) {
                output << tableau[0][num_rest + i] << " ";
            }
            output << endl;
        } else {
            output << "ilimitado" << endl;
            // Write a feasible solution
            //output << fixed << setprecision(DEC_PRECISION);
            int basic_row_index = -1;
            for (int i = 0; i < num_vars; ++i) {
                basic_row_index = -1;
                for (int j = 0; j < num_rest; ++j) {
                    if (basis[j] == 2 * num_rest + i) {
                        basic_row_index = j;
                        break;
                    }
                }
                if (basic_row_index != -1) {
                    output << tableau[basic_row_index + 1][total_cols - 1] << " ";
                } else {
                    output << "0.000 ";
                }
            }
            output << endl;
            // Write certificate of unboundedness
            for (int i = 0; i < num_vars; ++i) {
                basic_row_index = -1;
                if (2 * num_rest + i == col_pivot) {
                    output << "1.000 ";
                    continue;
                }
                for (int j = 0; j < num_rest; ++j) {
                    if (basis[j] == 2 * num_rest + i) {
                        basic_row_index = j;
                        break;
                    }
                }
                if (basic_row_index != -1) {
                    output << -tableau[basic_row_index + 1][col_pivot] << " ";
                } else {
                    output << "0.000 ";
                }
            }
            output << endl;
        }   
    }

    // Deallocate tableau
    for (int i = 0; i < num_rest + 2; ++i) {
        delete[] tableau[i];
    }
    delete[] tableau;

    // Deallocate basis
    delete[] basis;

    return 0;
}