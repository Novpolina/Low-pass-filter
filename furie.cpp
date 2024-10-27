
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <cassert> 
#include <fstream>

using namespace std;

const double all_time = 3e-4;
const double period = 1e-4;
const double T = 1e-5;
const int xm = 1;
//const double f0 = 40000;
const double A = 0.5;
const double dt = 0.000001;
const int diskret_num = all_time / dt;
double df = 0.5 / (1000 * T);



double K(double f, double f0) {
    return 1/sqrt(1+(1.0/f0)*(1.0/f0)*f*f);
}

double K1(double f, double f0) {
    if(f < f0){
        return 1;
    }
    else{
        return 0;
    }
}


// Function to calculate input signal
std::vector<double> signal(const std::vector<double>& t) {
    std::vector<double> ret(t.size());
    for (size_t i = 0; i < t.size(); ++i) {
       ret[i] = xm*std::fmod(t[i], period)/period + A*sin(8*M_PI*t[i]/T);
       // ret[i] = A*sin(8*M_PI*t[i]/T);

    }
    return ret;
}

// Функция для выполнения дискретного преобразования Фурье с заданными частотами
vector<complex<double>> dft(const vector<double>& input, const vector<double>& frequencies, const vector<double>& t) {
    int N = frequencies.size();  // Количество частот
    vector<complex<double>> output(N); // Вектор для хранения выходных данных

    for (int k = 0; k < N; k ++) {
        complex<double> sum(0, 0);
        for (size_t n = 0; n < input.size(); ++n) {
            sum += input[n] * std::exp(std::complex<double>(0, -2* M_PI * frequencies[k] * t[n])) * dt; // Суммирование
        }
        output[k] = sum;
    }
    return output;
}

// Функция для выполнения обратного дискретного преобразования Фурье с заданными частотами
vector<double> idft(const vector<complex<double>>& input, const vector<double>& frequencies, const vector<double>& t) {
    int N = frequencies.size();  // Количество частот
    vector<double> output(diskret_num); // Вектор для хранения выходных данных

    for (size_t n = 0; n < output.size(); ++n) {
        complex<double> sum(0, 0);
        for (int k = 0; k < static_cast<int>(input.size()) ; k ++ ) {
             sum += input[k] * std::exp(std::complex<double>(0, 2* M_PI * frequencies[k] * t[n])) * df; // Суммирование
        }
        complex<double> sum1 = sum * 1.0/( M_PI);
        output[n] = real(sum); // Нормализация
    }
    return output;
}

void save_to_file(const std::vector<double>& x, const std::vector<double>& y, const std::string& filename) {
    // Проверка, что оба вектора одинаковой длины
    if (x.size() != y.size()) {
        std::cerr << "Error: Vectors x and y must have the same length." << std::endl;
        return;
    }

    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Запись данных в файл
    for (size_t i = 0; i < x.size(); ++i) {
        outfile << x[i] << " " << y[i] << "\n";
    }

    outfile.close();
    std::cout << "Data successfully saved to " << filename << std::endl;
}

// Функция для проверки равенства массивов с заданной точностью
bool arraysEqual(const vector<double>& a, const vector<double>& b, double epsilon = 1e-10) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (abs(a[i] - b[i]) > epsilon) return false;
    }
    return true;
}
double GetFrequency(){
    double R, C;
    cout << "Enter value of R in Om \n";
    cin >> R;
    cout << "Enter value of C in F \n";
    cin >> C;
    double f0 = 1.0 /(2.0 * M_PI * R * C );
    return f0;
}

int main() {
    double f0 = GetFrequency();

    // Исходный сигнал
    std::vector<double> t(diskret_num);
    for (int i = 0; i < diskret_num; ++i) {
        t[i] = i * dt;
    }

    std::vector<double> input = signal(t);
    save_to_file(t, input, "x(t).txt" );


    // Заданное пространство частот
    std::vector<double> frequencies(10000);
    for (int i = 0; i < 10000; ++i) {
        frequencies[i] = i * df;
     }

    // Выполнение ДПФ
    vector<complex<double>> freqDomain = dft(input, frequencies, t);


    // Спектр входного сигнала
    std::vector<double> freqDomainReal(10000);

    for(int i = 0; i < freqDomain.size(); i++){
        freqDomainReal[i] = freqDomain[i].real(); 
    }
    save_to_file(frequencies, freqDomainReal, "spectrum.txt");


    // Преобразование сигнала ФНЧ

    for(int i = 0; i < frequencies.size(); i++){
        freqDomain[i] = freqDomain[i] * K(frequencies[i], f0);
    }

    //Спектр выхордного сигнала

    for(int i = 0; i < freqDomain.size(); i++){
        freqDomainReal[i] = freqDomain[i].real(); 
    }
      save_to_file(frequencies, freqDomainReal, "spectrum_filtered.txt");


    // Выполнение обратного ДПФ
    vector<double> reconstructedSignal = idft(freqDomain, frequencies, t);
    save_to_file(t, reconstructedSignal, "y(t).txt");

    // Проверка, совпадает ли восстановленный сигнал с исходным
    if (arraysEqual(input, reconstructedSignal)) {
        cout << "Signals are equal" << endl;
    } else {
        cout << "Signals are not equal" << endl;
    }

    return 0;
}


