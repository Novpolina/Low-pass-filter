import numpy as np
import matplotlib.pyplot as plt

# Чтение данных из файла
data = np.loadtxt('spectrum.txt')

# Извлечение значений частот и модулей преобразования Фурье
frequencies = data[:, 0]
magnitude = data[:, 1]

# Создание графика
plt.plot(frequencies, magnitude, label = "Spectrum", color='green')


# Чтение данных из файла
data = np.loadtxt('spectrum_filtered.txt')

# Извлечение значений частот и модулей преобразования Фурье
frequencies = data[:, 0]
magnitude = data[:, 1]

# Создание графика
plt.plot(frequencies, magnitude, label = "Spectrum_filtered", color='red')


plt.title('Спектр сигнала x(t)')
plt.xlabel('Частота (Гц)')
plt.ylabel('Сигнал (В)')
plt.grid(True)
plt.legend()
plt.show()