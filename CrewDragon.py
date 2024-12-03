import numpy as np
import matplotlib.pyplot as plt

g = 1.62     # Ускорение свободного падения на Луне
dt = 0.001   # Рассматриваем очень малые промежутки времени

class Apollo11:
	def __init__(self, M, fuelMass, V_p, consumption, H0, V0y, V_end):
		self.ApolloMass = M + fuelMass
		self.FuelMass = fuelMass
		self.FuelSpeed = V_p
		self.Consumption = consumption
		self.StartingHeight = H0
		self.StartingVelocity = V0y
		self.EndVelocity = V_end

	def HeightGraphic(self, _0_t0, _t0_t, H1, H2):
		fig, ax = plt.subplots(figsize=(8, 10))
		ax.plot(_0_t0, H1, color='black', label='Высота (без двигателя)')
		ax.plot(_t0_t, H2, color='orange', label='Высота (с двигателем)')
		ax.axhline(0, color='black', lw=0.5)
		ax.set_xlabel('Время (с)')
		ax.set_ylabel('Высота (м)')
		ax.set_title('Зависимость высоты от времени')
		ax.legend()
		ax.grid(True)
		plt.show()

	def AccelerationGraphic(self, _0_t0, _t0_t, a1, a2):
		fig, ax = plt.subplots(figsize=(8, 10))
		ax.plot(_0_t0, a1, color='black', label='Ускорение (без двигателя)')
		ax.plot(_t0_t, a2, color='orange', label='Ускорение (с двигателем)')
		ax.axhline(0, color='black', lw=0.5)
		ax.set_xlabel('Время (с)')
		ax.set_ylabel('Ускорение (м/с²)')
		ax.set_title('Зависимость ускорения от времени')
		ax.legend()
		ax.grid(True)
		plt.show()

	def VelocityGraphic(self, _0_t0, _t0_t, V1, V2):
		fig, ax = plt.subplots(figsize=(8, 10))
		ax.plot(_0_t0, V1, color='black', label='Скорость (без двигателя)')
		ax.plot(_t0_t, V2, color='orange', label='Скорость (с двигателем)')
		ax.set_xlabel('Время (с)')
		ax.set_ylabel('Скорость (м/с)')
		ax.set_title('Зависимость вертикальной скорости от времени')
		ax.legend()
		ax.grid(True)
		plt.show()

	def TimeCalculation(self):
		t0 = 0   # Время запуска двигателя
		t = 0    # Время посадки
		found = False

		# Находим время, при котором H = 0, решая квадратное уравнение
		t_end = (-self.StartingVelocity + (self.StartingVelocity ** 2 + 2 * g * self.StartingHeight) ** 0.5) / g

		# Поиск времени запуска двигателя
		for t0 in np.arange(0, t_end, dt):
			# Высота
			H1 = self.StartingHeight - self.StartingVelocity * t0 - 0.5 * g * t0 ** 2

			# Отрицательная высота не нужна
			if H1 < 0:
				break

			# Скорость
			V1 = self.StartingVelocity + g * t0
			t_max = (self.EndVelocity - V1) / (g - self.Consumption * self.FuelSpeed / self.ApolloMass)
			
			# Найдем t
			delta_t = t_max
			H2 = (H1 - V1 * delta_t - g * delta_t ** 2 / 2 + self.FuelSpeed * delta_t + self.FuelSpeed * self.ApolloMass / self.Consumption * (1 - self.Consumption * delta_t / self.ApolloMass) * np.log(1 - self.Consumption * delta_t / self.ApolloMass))
			
			if(H2 > 0):
				continue

			for t in np.arange(t0, t0 + t_max, dt):
				delta_t = t - t0
				remainingMass = self.FuelMass - self.Consumption * delta_t
				if remainingMass <= 0:
					break

				if delta_t > 0:
					H2 = (H1 - V1 * delta_t - g * delta_t ** 2 / 2 + self.FuelSpeed * delta_t + self.FuelSpeed * self.ApolloMass / self.Consumption * (1 - self.Consumption * delta_t / self.ApolloMass) * np.log(1 - self.Consumption * delta_t / self.ApolloMass))
					V2 = V1 + g * delta_t - self.FuelSpeed * np.log(self.ApolloMass / (self.ApolloMass - self.Consumption * delta_t))

					if H2 <= 0 and abs(V2) <= self.EndVelocity:
						found = True
						break
			if found:
				break

		if found:
			return t0, t
		else:
			print("Не удалось просчитать t0 и t. Возвращайтесь позже")


	def Prilunennie(self):
		t0, t = self.TimeCalculation()

		# Время от 0 до включения двигателя
		_0_t0 = np.linspace(0, t0, 100)
		# Время от включения двигателя до посадки
		_t0_t = np.linspace(t0, t, 100)

		# Ускорение до включения двигателя
		a1 = np.full_like(_0_t0, g)
		# Скорость от начала до включения двигателя
		V1 = self.StartingVelocity + g * _0_t0
		# Высота (путь) от начала до времени включения двигателя
		H1 = self.StartingHeight - self.StartingVelocity * _0_t0 - 0.5 * g * _0_t0 ** 2

		# Cкорость во время включения двигателя
		V = self.StartingVelocity + g * t0

		# Ускорение после включения двигателя
		a2 = g - (self.Consumption * self.FuelSpeed) / (self.ApolloMass - self.Consumption * (_t0_t - t0))
		# Скорость после включения двигателя
		V2 = (V - self.FuelSpeed * np.log(self.ApolloMass / (self.ApolloMass - self.Consumption * (_t0_t - t0))) + g * (_t0_t - t0))
		# Высота (путь)  после включения двигателя
		H2 = (H1[-1] - V1[-1] * (_t0_t - t0) - g * (_t0_t - t0) ** 2 / 2 + self.FuelSpeed * (_t0_t - t0) + self.FuelSpeed * self.ApolloMass / self.Consumption * (1 - self.Consumption * (_t0_t - t0) / self.ApolloMass) * np.log(1 - self.Consumption * (_t0_t - t0) / self.ApolloMass))

		print(f"Вертикальная скорость при h = 0 равна: {V2[-1]}")
		print(f"Необходимая высота для включения двигателя равна {H1[-1]}")

		self.HeightGraphic(_0_t0, _t0_t, H1, H2)
		self.AccelerationGraphic(_0_t0, _t0_t, a1, a2)
		self.VelocityGraphic(_0_t0, _t0_t, V1, V2)


starship = Apollo11(2150, 150, 3660, 15, 950, 61, 3)
starship.Prilunennie()




