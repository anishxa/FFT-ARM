# FFT-ARM
## FFT Application Using ARM Processor

#### Contributor: https://github.com/JoshLobo-21

We apply the FFT using a low-powered ARM core so that it is compact to carry and can be used in the field without much difficulty. Some notable applications of our project are:
1. Used for spying and eavesdropping (military applications)
2. Education and practical demonstrations
3. Signal repeaters and amplifiers

### Methodology & Observation:
Aimed to analyze the frequency components of the supplied analog signal. To achieve this, we utilized the STM32F401RE Nucleo board and C code. The analog signal was collected using the onboard ADC with a potentiometer as a substitute. The data from the ADC was then transmitted via UART to the console and fed into our FFT code, which computed and displayed the magnitude and frequency components of the waveform. By employing a radix-4 butterfly and 2 stages, we efficiently computed the FFT of a 16-point sample, unlike the radix-2 FFT, which requires 4 stages for the same task.

### Result:
We successfully computed the FFT of a 16-point sample and validated the results by comparing them with MATLAB's built-in function. Additionally, we plotted the waveform and frequency amplitude graph using MATLAB. Since it is at a small scale, it can be computed further and applied for audio processing, vibration analysis and so on.
