#include <iostream>
#include <omp.h>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>

#pragma GCC optimize("O3")

using namespace std;

int main(int argc, char *argv[]) {
    double threshold = atof(argv[4]);
    int threads_num = atoi(argv[1]);
#ifdef _OPENMP
    omp_set_num_threads(threads_num);
#endif
    ifstream in(argv[2], ios::binary);
    ofstream out(argv[3]);
    in.tie(0);
    out.tie(0);
    if (!in.is_open()) {
        cerr << "File Not Found" << endl;
        exit(1);
    }
    int k = 1;
    string cur="";
    string header[4];
    getline(in, header[0]);
    string line;
    getline(in, line);
    for (int i = 0; i < line.length(); i++) {
        if (isspace(line[i])) {
            header[1] = line.substr(0, i);
            header[2] = line.substr(i+1, line.length() - i - 1);
            break;
        }
    }
    getline(in, header[3]);
    if(isspace(in.peek())) {
        in.get();
    }
    cout << header[0] << endl << header[1] << " " << header[2] << endl <<header[3]<<endl;
    if (header[0] != "P5" && header[0] != "P6" || header[3] != "255") {
        cerr << "Invalid type" << endl;
        exit(2);
    }
    if (header[0] == "P6") {
        k = 3;
    }
    out << header[0] << endl << header[1] << " " << header[2] << endl << 255<<endl;
    int *buffer;
    int i = 0, size = k * atoi(header[1].c_str()) * atoi(header[2].c_str());
    buffer = new int[size];
    while (i < size) {
        buffer[i]=in.get();
        i++;
    }
    auto start = chrono::high_resolution_clock::now();
    int cnt_pixels[3][256] = {0};
    int maxes [3]{0,0,0};
    int minis [3]{255,255,255};
#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(none) shared(k,size, buffer, threads_num, cnt_pixels)
#endif
    for (int ch = 0; ch < k; ch++) {
        for (int i = ch; i < size; i += k) {
            cnt_pixels[ch][buffer[i]]++;
        }
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(none) shared(k,size, minis, maxes, buffer, threshold, threads_num, cnt_pixels)
#endif
    for (int ch = 0; ch < k; ch++) {
        int left_sum = 0, right_sum = 0;
        bool flag_min = true, flag_max = true;
        for (int i = 0; i < 256; i++) {
            left_sum += cnt_pixels[ch][i];
            if (left_sum >= size / k * threshold && flag_min) {
                minis[ch] = i;
                flag_min = false;
            }
            right_sum += cnt_pixels[ch][255 - i];
            if (right_sum >= size / k * threshold && flag_max) {
                maxes[ch] = 255 - i;
                flag_max = false;
            }
        }
    }
    for (int ch = 0; ch < k; ch++) {
        if (maxes[ch] < minis[ch]) {
            maxes[ch] = (minis[ch] + maxes[ch]) / 2;
            minis[ch] = maxes[ch];
        }
    }
    int *final = new int[256];
    int minimum = min(minis[0], min(minis[1], minis[2]));
    int maximum = max(maxes[0], max(maxes[1], maxes[2]));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(none) shared(final, maximum, minimum)
#endif
    for (int i = 0; i < 256; i++) {
        if (minimum == maximum) {
            final[i] = minimum;
        } else if (i <= minimum) {
            final[i] = 0;
        } else if (i >= maximum) {
            final[i] = 255;
        } else {
            final[i] = (int) (255.0 * ((double) i - minimum) / (maximum - minimum));
        }
    }
    auto stop = chrono::high_resolution_clock::now();
    auto total_time = chrono::duration_cast<chrono::milliseconds>(stop - start).count();
    int pos = 0;
    while (pos < size) {
        out << (unsigned char) final[buffer[pos]];
        pos++;
    }
    out.close();
    delete [] final;
    delete [] buffer;
    fprintf(stdout, "Time (%i thread(s)): %i ms\n", threads_num, total_time);
    return 0;
}
