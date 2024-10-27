import numpy as np
import matplotlib.pyplot as plt

# Чтение данных из файла
data = np.loadtxt('x(t).txt')

# Извлечение значений частот и модулей преобразования Фурье
time_x = data[:, 0]
signal_x = data[:, 1]

# Создание графика
plt.plot(time_x, signal_x, label='Входной сигнал', color='green')

# Чтение данных из файла
data = np.loadtxt('y(t).txt')

# Извлечение значений частот и модулей преобразования Фурье
time_y = data[:, 0]
signal_y = data[:, 1]

# Создание графика
plt.plot(time_y, signal_y, label = "Выходной сигнал", color='red')


plt.title('Преобразование сигнала фильтром низких частот')
plt.xlabel('Время')
plt.ylabel('Сигнал')
plt.grid(True)
plt.legend()
plt.show()
