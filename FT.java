package priv.droptea.emotion.ft;

import priv.droptea.emotion.util.ChartTool;

public class FT {
    public static void main(String[] args) {
    	//模拟一个由两个频率不同的正弦波组成的音频数据
    	int N = 5000;
    	//5个点一个周期，假设这N个数据是Ts秒的音频数据，那么f=5000/5*1/Ts=1000/Ts,由k=f*Ts,所以傅里叶变换后k=1000的值就是这个sin函数的振幅
    	double f_sin1 = 1.0/5;
    	//改变第一个正弦函数的相位
    	int phase_sin1 = 1;
    	//振幅
    	double amplitude_sin1 = 0.8;
    	//同上，傅里叶变换后k=500的值就是这个sin函数的振幅
    	double f_sin2 = 1.0/10;
    	//振幅
    	double amplitude_sin2 = 0.4;
    	Complex[] data = Complex.initArray(N);
    	for (int i = 0; i < data.length; i++) {
    		data[i].re=amplitude_sin1*Math.sin(2*Math.PI*i*f_sin1+phase_sin1)+amplitude_sin2*Math.sin(2*Math.PI*i*f_sin2);//0.8*Math.sin(2*PI*k1*i);
		}
    	new ChartTool("原函数",modulus(data));
    	//dftForward
    	long timeDftForward = System.currentTimeMillis();
    	Complex[] resultDftForward = dftForward(data);
    	new ChartTool("dftForward",modulus(resultDftForward));
    	System.out.println("runtime_dftForward:"+(System.currentTimeMillis()-timeDftForward));
    	//dftInverse
    	long timeDftInverse = System.currentTimeMillis();
    	Complex[] resultDftInverse =  dftInverse(resultDftForward);
    	new ChartTool("dftInverse",modulus(resultDftInverse));
    	System.out.println("runtime_dftInverse:"+(System.currentTimeMillis()-timeDftInverse));
    	//fftForward
    	long timeFftForward = System.currentTimeMillis();
    	Complex[] resultFftForward = fftForward(data);
    	new ChartTool("fftForward",modulus(resultFftForward));
    	System.out.println("runtime_fftForward:"+(System.currentTimeMillis()-timeFftForward));
    	//fftForward
    	long timeFftInverse = System.currentTimeMillis();
    	Complex[] resultFftInverse = fftInverse(resultFftForward);
    	new ChartTool("fftInverse",modulus(resultFftInverse));
    	System.out.println("runtime_fftInverse:"+(System.currentTimeMillis()-timeFftInverse));
    	
    	new ChartTool("相位谱",atan(resultFftForward));
	}
    
    public static Complex[] fftForward(Complex[] x) {
		int N = x.length;
		// 因为exp(-2i*n*PI)=1，n=1时递归原点
		if (N == 1){
			return x;
		}
		// 如果信号数为奇数，使用dft计算
		if (N % 2 != 0) {
			return dftForward(x);
		}
		// 提取下标为偶数的原始信号值进行递归fft计算
		Complex[] even = new Complex[N / 2];
		for (int k = 0; k < N / 2; k++) {
			even[k] = x[2 * k];
		}
		Complex[] evenValue = fftForward(even);
		// 提取下标为奇数的原始信号值进行fft计算
		// 节约内存
		Complex[] odd = even;
		for (int k = 0; k < N / 2; k++) {
			odd[k] = x[2 * k + 1];
		}
		Complex[] oddValue = fftForward(odd);
		// 偶数+奇数
		Complex[] result = new Complex[N];
		Complex temp = new Complex();
		for (int k = 0; k < N / 2; k++) {
			double p = -2*Math.PI*k/N;
			temp.re=Math.cos(p);
			temp.im=Math.sin(p);
			Complex multValue = complexMult(temp,oddValue[k]);
			result[k] =complexAdd(evenValue[k],multValue);
			result[k + N/2] = complexMinus(evenValue[k],multValue);
		}
		return result;
	}
    
    public static Complex[] fftInverse(Complex[] x) {
		int N = x.length;
		// 因为exp(-2i*n*PI)=1，n=1时递归原点
		if (N == 1){
			return x;
		}
		// 如果信号数为奇数，使用dft计算
		if (N % 2 != 0) {
			return dftForward(x);
		}
		// 提取下标为偶数的原始信号值进行递归fft计算
		Complex[] even = new Complex[N / 2];
		for (int k = 0; k < N / 2; k++) {
			even[k] = x[2 * k];
		}
		Complex[] evenValue = fftForward(even);
		// 提取下标为奇数的原始信号值进行fft计算
		// 节约内存
		Complex[] odd = even;
		for (int k = 0; k < N / 2; k++) {
			odd[k] = x[2 * k + 1];
		}
		Complex[] oddValue = fftForward(odd);
		// 偶数+奇数
		Complex[] result = new Complex[N];
		Complex temp = new Complex();
		for (int k = 0; k < N / 2; k++) {
			double p = 2*Math.PI*k/N;
			temp.re=Math.cos(p);
			temp.im=Math.sin(p);
			Complex multValue = complexMult(temp,oddValue[k]);
			result[k] =complexAdd(evenValue[k],multValue);
			result[k].re /= N;
			result[k].im /= N;
			result[k + N/2] = complexMinus(evenValue[k],multValue);
			result[k + N/2].re /= N;
			result[k + N/2].im /= N;
		}
		return result;
	}

	//离散傅里叶变换
	public static Complex[] dftForward(Complex[] input) {
		Complex[] output =  Complex.initArray(input.length);
		Complex temp = new Complex();
		int N = input.length;
		for(int k=0;k<N;k++) {
			for(int n=0;n<N;n++) {
				double p = -2*Math.PI/N*n*k;
				temp.re=Math.cos(p);
				temp.im=Math.sin(p);
				output[k] = complexAdd(output[k],complexMult(input[n],temp));
			}
		}
		return output;
	}
	//离散傅里叶逆变换
	public static Complex[] dftInverse(Complex[] input) {
		Complex[] output =  Complex.initArray(input.length);
		Complex temp = new Complex();
		int N = input.length;
		for(int k=0;k<N;k++) {
			for(int n=0;n<N;n++) {
				temp.re=Math.cos(2*Math.PI/N*n*k);
				temp.im=Math.sin(2*Math.PI/N*n*k);
				output[k] = complexAdd(output[k],complexMult(input[n],temp));
			}
			output[k].re /= N;
			output[k].im /= N;
		}
		return output;
	}
	
	public static double[] modulus(Complex[] complaxs) {
		double[] result = new double[complaxs.length];
		for (int i = 0; i <complaxs.length ; i++) {
			result[i] = Math.sqrt(complaxs[i].re*complaxs[i].re+complaxs[i].im*complaxs[i].im);
		}
		return result;
	}
	public static double[] atan(Complex[] complaxs) {
		double[] result = new double[complaxs.length];
		for (int i = 0; i <complaxs.length ; i++) {
			if(Math.sqrt(complaxs[i].re*complaxs[i].re+complaxs[i].im*complaxs[i].im)>0.2) {
				result[i] = Math.atan(complaxs[i].im/complaxs[i].re);
			}
		}
		return result;
	}
	
	public static Complex complexAdd(Complex a, Complex b){ //复数加
		Complex rt = new Complex();
	    rt.re = a.re + b.re;
	    rt.im = a.im + b.im;
	    return rt;
	}
	public static Complex complexMinus(Complex a, Complex b){ //复数加
		Complex rt = new Complex();
	    rt.re = a.re - b.re;
	    rt.im = a.im - b.im;
	    return rt;
	}
	 
	public static Complex complexMult(Complex a, Complex b){ //复数乘
		Complex rt= new Complex();
	    rt.re = a.re*b.re-a.im*b.im;
	    rt.im = a.im*b.re+a.re*b.im;
	    return rt;
	}
	public static class Complex{
		public double re;
		public double im;
		public Complex() {
		}
		public Complex(double re,double im) {
			this.re = re;
			this.im = im;
		}
		public static Complex[] initArray(int N) {
			Complex[] complaxs = new Complex[N];
			for (int i = 0; i < N; i++) {
				complaxs[i] = new Complex();
			}
			return complaxs;
		}
	}
}
