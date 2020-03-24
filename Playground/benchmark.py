from numpy import *
import time

signal = random.random_sample(4096)

prewarm = 3
loop = 100000

for i in range(0, prewarm):
    fft.fft(signal)

now = time.time()
for i in range(0, loop):
    fft.fft(signal)
print(str((time.time() - now) / loop) + " seconds")
