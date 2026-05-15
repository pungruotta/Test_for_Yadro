#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include <iostream>

using namespace std;

const double PI = acos(-1);

class FastFFT {
public:
    static int Power(int n){
        int k_2 = n;
        while (k_2%2==0){
            k_2 = k_2/2;
        }
        if (k_2 == 1)
            return 2;

        int k_3 = n;
        while (k_3%3==0){
            k_3 = k_3/3;
        }
        if (k_3 == 1)
            return 3;
        
        int k_5 = n;
        while (k_5%5==0){
            k_5 = k_5/5;
        }
        if (k_5 == 1)
            return 5;
        
    }
 
    static void fft (vector<complex<double>> & a, bool invert) {
    	int n = (int) a.size();
	    if (n == 1)  return;
        int b = Power(n);
            if(b==2){
 
    	        vector<complex<double>> a0 (n/2),  a1 (n/2);
    	        for (int i=0, j=0; i<n; i+=2, ++j) {
    		        a0[j] = a[i];
	    	        a1[j] = a[i+1];
    	        }
	            fft (a0, invert);
    	        fft (a1, invert);
 
    	        double ang = 2*PI/n * (invert ? -1 : 1);
    	        complex<double> w (1),  wn (cos(ang), sin(ang));
	            for (int i=0; i<n/2; ++i) {
	    	        a[i] = a0[i] + w * a1[i];
	    	        a[i+n/2] = a0[i] - w * a1[i];
	    	        if (invert)
	    		        a[i] /= 2,  a[i+n/2] /= 2;
	    	        w *= wn;
        
      	        }
            }else if(b==3){
                vector<complex<double>> a0(n / 3), a1(n / 3), a2(n / 3);
                for (int i = 0; 3 * i < n; i++) {
                    a0[i] = a[3 * i];
                    a1[i] = a[3 * i + 1];
                    a2[i] = a[3 * i + 2];
                }
                fft(a0, invert);
                fft(a1, invert);
                fft(a2, invert);

                double angle = 2 * PI / n * (invert ? -1 : 1);
                complex<double> w_n = complex<double>(cos(angle), -sin(angle));
                complex<double> w = 1;

                for (int i = 0; i < n / 3; i++) {
                    complex<double> twid1 = w;
                    complex<double> twid2 = w * w;

                    complex<double> p0 = a0[i];
                    complex<double> p1 = twid1 * a1[i];
                    complex<double> p2 = twid2 * a2[i];
                    double ang = 2 * PI / 3 * (invert ? -1 : 1);
                    complex<double> w3_1 = complex<double>(cos(ang), -sin(ang));
                    complex<double> w3_2 = complex<double>(cos(2 * ang), -sin(2 * ang));

                    a[i] = p0 + p1 + p2;
                    a[i + n / 3] = p0 + p1 * w3_1 + p2 * w3_2;
                    a[i + 2 * n / 3] = p0 + p1 * w3_2 + p2 * w3_1;
                    if (invert)
                        a[i] /= 3,  a[i+n/3] /= 3, a[i + 2 * n / 3] /= 3;

                    w *= w_n;
                }

                
            }else if(b==5){
                vector<complex<double>> a0(n / 5), a1(n / 5), a2(n / 5), a3(n / 5), a4(n / 5);
                for (int i = 0; 5 * i < n; i++) {
                    a0[i] = a[5 * i];
                    a1[i] = a[5 * i + 1];
                    a2[i] = a[5 * i + 2];
                    a3[i] = a[5 * i + 3];
                    a4[i] = a[5 * i + 4];
                }
                fft(a0, invert);
                fft(a1, invert);
                fft(a2, invert);
                fft(a3, invert);
                fft(a4, invert);
            
                double angle = 2 * PI / n * (invert ? -1 : 1);
                complex<double> w_n = complex<double>(cos(angle), -sin(angle));
                complex<double> w = 1;

                double ang = 2 * PI / 5 * (invert ? -1 : 1);
                complex<double> w5_1 = complex<double>(cos(ang), -sin(ang));
                complex<double> w5_2 = complex<double>(cos(2 * ang), -sin(2 * ang));
                complex<double> w5_3 = complex<double>(cos(3 * ang), -sin(3 * ang));
                complex<double> w5_4 = complex<double>(cos(4 * ang), -sin(4 * ang));
                for (int i = 0; i < n / 5; i++) {
                    complex<double> twid1 = w;
                    complex<double> twid2 = twid1 * w;
                    complex<double> twid3 = twid2 * w;
                    complex<double> twid4 = twid3 * w;

                    complex<double> p0 = a0[i];
                    complex<double> p1 = twid1 * a1[i];
                    complex<double> p2 = twid2 * a2[i];
                    complex<double> p3 = twid3 * a3[i];
                    complex<double> p4 = twid4 * a4[i];

                    a[i] = p0 + p1 + p2 + p3 + p4;
                    a[i + n / 5] = p0 + p1 * w5_1 + p2 * w5_2 + p3 * w5_3 + p4 * w5_4;
                    a[i + 2 * n / 5] = p0 + p1 * w5_2 + p2 * w5_4 + p3 * w5_1 + p4 * w5_3;
                    a[i + 3 * n / 5] = p0 + p1 * w5_3 + p2 * w5_1 + p3 * w5_4 + p4 * w5_2;
                    a[i + 4 * n / 5] = p0 + p1 * w5_4 + p2 * w5_3 + p3 * w5_2 + p4 * w5_1;

                    if (invert) 
                        a[i] /= 5,  a[i+n/5] /= 5, a[i + 2 * n / 5] /= 5, a[i + 3 * n / 5] /= 5, a[i + 4 * n / 5] /= 5;
                    w *= w_n;
                }

                
            }
    
    
    }
};


int main(){
    int n = 8;
    
    vector<complex<double>> data(n);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(-1.0, 1.0);
    for (auto& c : data) {
        c = complex<double>(dis(gen), dis(gen));
    }

    vector<complex<double>> original_data = data;

    FastFFT::fft(data, false);

    cout << "\nTransformed data:\n";
    for (const auto& c : data)
        cout << c << "\n";

    FastFFT::fft(data, true);

    double error_sum = 0.0;
    for (size_t i = 0; i < n; ++i) {
        error_sum += norm(original_data[i] - data[i]);
    }
    double mse = error_sum / n;

    cout << "Original data:\n";
    for (const auto& c : original_data)
        cout << c << "\n";

    cout << "\nTransformed back data:\n";
    for (const auto& c : data)
        cout << c << "\n";
    cout << "\nError sum:\n";
    cout << mse;
    return 0;
}