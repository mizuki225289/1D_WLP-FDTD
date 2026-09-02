/* ********************
2次元・3次元の連続したメモリ領域を確保・開放するための関数
*********************** */

double** allocate_memory2d(const int M, const int N, const double init_value) {
    double *v1d = new double [M*N]; /*1次元のまとまった領域確保　先頭のアドレスがv1d*/
    double **v = new double* [M]; /*M*N 配列の先頭のアドレスを入れておく配列（要素数M）を確保 先頭のアドレスがv*/
    
    for(int m=0; m < M; m++) {
        v[m] = v1d + m * N; /*M*N 配列の先頭のアドレスを入れておく配列（要素数M）に各行のアドレスを入れる*/
        for(int n=0; n < N;n++) {
            v[m][n] = init_value;
        }
    }

    return v;
}

void deallocate_memory2d(double** v) {
    delete [] v[0];
    delete [] v;
}

double*** allocate_memory3d(const int P, const int Q, const int R, double init_value) {
    double ***v = new double** [P];

    for(int p=0; p < P; p++) {
        v[p] = new double*[Q];
        for(int q=0; q < Q; q++) {
            v[p][q] = new double[R];
            for(int r=0; r < R; r++) {
                v[p][q][r] = init_value;
            }
        }
    }

    return v;
}

void deallocate_memory3d(double*** v, const int P, const int Q) {
    if (v == nullptr) {
        return;
    }

    for (int p = 0; p < P; ++p) {
        for (int q = 0; q < Q; ++q) {
            delete[] v[p][q];
        }
        delete[] v[p];
    }
    delete[] v;
}