from scipy.io import wavfile
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.signal import spectrogram
from scipy.signal import lfilter
from scipy.signal import tf2zpk
from scipy.signal import freqz

def fisrt(fs, data):
    
    print(len(data))#pocet vzorkov
    print(data.min())#min hodnota
    print(data.max())#max hodnota
    print(data.size/fs)#casova dlzka signalu na sekundy
    
    #vykreslenie signalu
    t = np.arange(data.size) / fs
    plt.figure(figsize=(6,3))
    plt.plot(t, data)
    plt.gca().set_xlabel('$t[s]$')
    plt.gca().set_title('Úloha 1')
    plt.tight_layout()
    plt.show()

def second(fs, data):
    #ustrednenie
    data = data-np.mean(data)#od dat odcitam ich priemernu hodnotu
    
    #normalizacia
    data_max = max(abs(data))
    data = data / data_max
    
    count = 0
    array = np.array([])
    matrix = np.array([])

    #vzorkovanie
    i = 0
    while i!=data.size:
        if count <=1023:
            array = np.insert(array,count,data[i])
            if count == 1023:
                matrix = np.append(matrix,array)
                array = np.array([])
                i -= 512
                count = -1
        i+=1
        count += 1               
    matrix = np.reshape(matrix, (-1, 1024))
       
    t = np.arange(matrix[11].size) / fs
    t = t + ((data.size/fs / 61)*11)#prepocitanie casu na grafe na jeho skutocny cas

    plt.figure(figsize=(6,3))
    plt.plot(t, matrix[11])
    plt.gca().set_xlabel('$t[s]$')
    word = "Úloha 2" 
    plt.gca().set_title(word)
    plt.tight_layout()
    plt.show()
    return matrix[11]#ramec pouzijem aj v ulohe 3

def third(fs, data):
    
    n = np.arange(data.size)
    k = n.reshape(data.size, 1)
    e = np.exp(-1j *2* np.pi * k * n / data.size)
    X = np.dot(e, data)#nasobenie bazi a matice dat
    X = np.abs(X)
    X = X[:512]#zobrazenie len polovice dft
    freq = np.arange(0, fs/2, fs/1024)
    plt.figure(figsize=(6,3))
    plt.plot(freq, X)
    plt.gca().set_xlabel('$frekvencia[Hz]$')
    word = "Moja DFT" 
    plt.gca().set_title(word)
    plt.tight_layout()
    plt.show()
    
    data = np.fft.fft(data)
    data = np.abs(data)
    data = data[:512]
    freq = np.arange(0, fs/2, fs/1024)
    plt.figure(figsize=(6,3))
    plt.plot(freq, data)
    plt.gca().set_xlabel('$frekvencia[Hz]$')
    word = "Funkcia np.fft.fft" 
    plt.gca().set_title(word)
    plt.tight_layout()
    plt.show()

    print(np.allclose(X, data))

def fourth(fs, data):
    #ustrednenie
    data = data-np.mean(data)
    #normalizacia
    data_max = max(abs(data))
    data = data / data_max
    matrix1 = np.array([])
    timeee = data.size/fs
    freq = fs/2
    
    trash,trash1,matrix1 = spectrogram(data,fs=8000,nperseg=1024,noverlap=512)
    matrix1 = 10 * np.log10(matrix1)

    plt.figure(figsize=(9,3))
    plt.imshow(matrix1, origin="lower",aspect="auto",extent=[0,timeee,0,freq])
    plt.gca().set_xlabel('Čas [s]')
    plt.gca().set_ylabel('Frekvencia [Hz]')
    plt.gca().set_title("Úloha 4")
    cbar = plt.colorbar()
    cbar.set_label('Spektrálna hustota výkonu [dB]', rotation=270, labelpad=15)

    plt.tight_layout()
    plt.show()


def fifth():

    r_val = [725, 1460, 2175, 2900]#hodnoty zistene manualne zo spektrogramu
    expected_val = [725, 725*2, 725*3, 725*4]#predpokladane hodnoty
    for i in range(len(r_val)):
        if r_val[i] == expected_val[i]:
            print("hodnota",i+1, "sa zhoduje s predpokladanou hodnotou")
        else:
            print("hodnota", i+1, "sa nezhoduje s predpokladanou hodnotou o", abs(expected_val[i] - r_val[i]))
 

def sixth(fs, data):

    t = np.array([])
    i = 0
    while i!= len(data):
        t = np.append(t,i/fs) #hodnoty ake ma cosinusovka v danom case
        i+=1

    #tvorenie cosinusoviek
    cos1 = np.cos(2 * np.pi * 725 * t)
    cos2 = np.cos(2 * np.pi * 1460 * t)
    cos3 = np.cos(2 * np.pi * 2175 * t)
    cos4 = np.cos(2 * np.pi * 2900 * t)
    
    output = np.array([])
    output = cos1+cos2+cos3+cos4#finalna cosinusovka
    #ustrednenie
    output = output-np.mean(output)#od dat odcitam ich priemernu hodnotu
    
    #normalizacia
    output_max = max(abs(output))
    output = output / output_max

    wavfile.write('../audio/ 4cos.wav',16000,(output * np.iinfo(np.int16).max).astype(np.int16))
    
    matrix1 = np.array([])
    tm = data.size/fs
    trash,trash1,matrix1 = spectrogram(output,fs=8000,nperseg=1024, noverlap=512)
    matrix1 = 10 * np.log10(matrix1)
    plt.figure(figsize=(9,3))
    plt.imshow(matrix1, extent=[0,tm,0,fs/2], origin="lower",aspect="auto")
    plt.gca().set_xlabel('Čas [s]')
    plt.gca().set_ylabel('Frekvencia [Hz]')
    plt.gca().set_title("Spektrogram")
    cbar = plt.colorbar()
    cbar.set_label('Spektrálna hustota výkonu [dB]', rotation=270, labelpad=15)

    plt.tight_layout()
    plt.show()

def seventh(fs):

    n = np.array([])
    frequency = [725,1460,2175,2900]
    #vypocet koeficientov filtra
    for i in range(1,5):
        w = 2*math.pi*frequency[i-1]/fs
        n = np.append(n,math.e**(1j*w))
    

    nkz = np.array([])
    nkz = np.conjugate(n)
    
    filter1 = np.array([])
    filter1 = np.append(filter1, n)
    filter1 = np.append(filter1, nkz)


    filter1 = np.poly(filter1)

    print("koeficienty su:",filter1)

    #pocet bodov v impulznej odozve
    N_imp = 9
    imp = [1, *np.zeros(N_imp-1)]
    h = lfilter(filter1, [1], imp)

    plt.figure(figsize=(5,3))
    plt.stem(np.arange(N_imp), h, basefmt=' ')
    plt.gca().set_xlabel('$n$')
    plt.gca().set_title('Impulzná odozva $h[n]$')

    plt.grid(alpha=0.5, linestyle='--')

    plt.tight_layout()
    plt.show()
    return filter1
   

def eighth(filter):
    z, p, k = tf2zpk(filter, [1])

    ang = np.linspace(0, 2*np.pi,100)
    plt.figure(figsize=(4,3.5))
    plt.plot(np.cos(ang), np.sin(ang))

    #zobrezenie nul a polov
    plt.scatter(np.real(z), np.imag(z), marker='o', facecolors='none', edgecolors='r', label='nuly')
    plt.scatter(np.real(p), np.imag(p), marker='x', color='g', label='póly')

    plt.gca().set_xlabel('Reálna zložka $\mathbb{R}\{$z$\}$')
    plt.gca().set_ylabel('Imaginárna zložka $\mathbb{I}\{$z$\}$')

    plt.grid(alpha=0.5, linestyle='--')
    plt.legend(loc='upper right')

    plt.tight_layout()
    plt.show()


def ninth(fs, filter):

    w, H = freqz(filter, [1])
    
    plt.figure(figsize=(6,3))
    plt.plot(w / 2 / np.pi * fs, np.abs(H))
    plt.gca().set_xlabel('Frekvencia [Hz]')
    plt.gca().set_title('Úloha 9 Modul frekvenčnej charakteristiky $|H(e^{j\omega})|$')

    plt.tight_layout()
    plt.show()

def tenth(fs, data, filter):
    #ustrednenie
    data = data-np.mean(data)#od dat odcitam ich priemernu hodnotu
    
    #normalizacia
    data_max = max(abs(data))
    data = data / data_max
    
    sf = lfilter(filter, [1], data)#aplikacia filtra na signal
    sf_max = max(abs(sf))
    sf = sf / sf_max #normalizacia
    
    t = np.arange(sf.size) / fs
    plt.figure(figsize=(6,3))
    plt.plot(t, sf)
    plt.gca().set_xlabel('$t[s]$')
    plt.gca().set_ylabel('$frekvencia[Hz]$')
    plt.gca().set_title('Úloha 10')
    plt.tight_layout()
    plt.show()

    wavfile.write('../audio/clean_z.wav',16000,(sf * np.iinfo(np.int16).max).astype(np.int16))




fs, data = wavfile.read('../audio/xpolon03.wav')#nacitanie vstupneho signalu

fisrt(fs, data)#uloha 1
matrix1 = second(fs, data)#uloha 2
third(fs,matrix1)#uloha 3
fourth(fs, data)#uloha 4
fifth()#uloha 5
sixth(fs, data)#uloha 6
filter1= seventh(fs)#uloha 7
eighth(filter1)#uloha 8
ninth(fs, filter1)#uloha 9
tenth(fs, data, filter1)#uloha 10




