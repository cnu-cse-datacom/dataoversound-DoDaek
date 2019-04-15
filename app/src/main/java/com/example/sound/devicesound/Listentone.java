package com.example.sound.devicesound;

import android.media.AudioFormat;
import android.media.AudioRecord;
import android.media.MediaRecorder;
import android.util.Log;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.*;

import java.util.ArrayList;
import java.util.List;

import google.zxing.common.reedsolomon.GenericGF;
import google.zxing.common.reedsolomon.ReedSolomonDecoder;
import google.zxing.common.reedsolomon.ReedSolomonException;

import static java.lang.Math.abs;
import static java.lang.Math.min;

public class Listentone {

    int HANDSHAKE_START_HZ = 4096;
    int HANDSHAKE_END_HZ = 5120 + 1024;

    int START_HZ = 1024;
    int STEP_HZ = 256;
    int BITS = 4;

    int FEC_BYTES = 4;

    private int mAudioSource = MediaRecorder.AudioSource.MIC;     // mic = alsaaudio.PCM(alsaaudio.PCM_CAPTURE, alsaaudio.PCM_NORMAL, device="default")
    private int mSampleRate = 44100;                              // mic.setrate(44100)
    private int mChannelCount = AudioFormat.CHANNEL_IN_MONO;      // mic.setchannels(1)
    private int mAudioFormat = AudioFormat.ENCODING_PCM_16BIT;    // mic.setformat(alsaaudio.PCM_FORMAT_S16_LE)
    private float interval = 0.1f;                                // interval=0.1

    private int mBufferSize = AudioRecord.getMinBufferSize(mSampleRate, mChannelCount, mAudioFormat);

    public AudioRecord mAudioRecord = null;
    int audioEncodig;
    boolean startFlag;
    FastFourierTransformer transform;


    public Listentone(){

        transform = new FastFourierTransformer(DftNormalization.STANDARD);                                      // STANDARD 옵션 : FFT의 스케일을 정규화하지 않겠다는 것
        startFlag = false;                                                                                      // in_packet = False
        mAudioRecord = new AudioRecord(mAudioSource, mSampleRate, mChannelCount, mAudioFormat, mBufferSize);
        mAudioRecord.startRecording();                                                                          // l, data = mic.read()

    }

    private int findPowerSize(int value){
        int num = 1;
        while(true){
            if(value < num){
                if(num - value > value - (num/2)){
                    return num/2;
                }
                return num;
            }
            num *= 2;
        }
    }

    private double findFrequency(double[] toTransform){
        int len = toTransform.length;                                                  // 매개변수인 배열의 길이
        double[] real = new double[len];                                               // 실수부를 담을 배열
        double[] img = new double[len];                                                // 허수부를 담을 배열
        double realNum;                                                                // 실수부
        double imgNum;                                                                 // 허수부
        double[] mag = new double[len];                                                // 크기를 담을 배열 : sin 함수에서 A의 값을 저장

        Complex[] complx = transform.transform(toTransform, TransformType.FORWARD);    // ==> w
        Double[] freq = this.fftfreq(complx.length, 1);                        // ==> freqs

        for(int i = 0 ; i < complx.length ; i++){
            realNum = complx[i].getReal();                                              // toTransform 배열에 대해 FFT 실행 후의 실수부
            imgNum = complx[i].getImaginary();                                          // toTransform 배열에 대해 FFT 실행 후의 허수부
            mag[i] = Math.sqrt((realNum * realNum) + (imgNum * imgNum));                // A의 값 구하기
        }

        int peak_coeff = 0;                   // peak frequency 의 index
        double peak = abs(mag[0]);            // peak frequency
        for(int i = 1 ; i < complx.length ; i++){
            if(abs(mag[i]) > peak){           // 현재 peak frequency 보다 크다면
                peak_coeff = i;               // peak frequency 의 index 변경
                peak = abs(mag[i]);           // peak frequency 저장
            }
        }

        double peakFreq = freq[peak_coeff];
        double result = peakFreq * mSampleRate;

        return abs(result);
    }

    private Double[] fftfreq(int length, int duration){
        double val = 1.0 / (length * duration);
        Double[] results = new Double[length];
        int[] tmp = new int[length];
        int N = (length - 1) / 2 + 1;
        for(int i = 0 ; i < N ; i++){
            tmp[i] = i;
        }
        for(int i = N, p2 = -(length / 2); i < length ; i++, p2++){
            tmp[i] = p2;
        }
        for(int i = 0 ; i < length ; i++){
            results[i] = tmp[i] * val;
        }
        return results;
    }
    /*
    def fftfreq(n, d=1.0):
        if not isinstance(n, integer_types):
            raise ValueError("n should be an integer")
        val = 1.0 / (n * d)
        results = empty(n, int)
        N = (n-1)//2 + 1
        p1 = arange(0, N, dtype=int)
        results[:N] = p1
        p2 = arange(-(n//2), 0, dtype=int)
        results[N:] = p2
        return results * val
     */

    public void PreRequest(){

        int blocksize = findPowerSize((int)(long)Math.round(interval/2*mSampleRate));
        short[] buffer = new short[blocksize];
        double[] transformed = new double[blocksize];
        List<Double> packet = new ArrayList<>();

        while(true){
            int bufferedReadResult = mAudioRecord.read(buffer, 0, blocksize);
            //Log.d("bufferedReadResult : ", Integer.toString(bufferedReadResult));
            if(bufferedReadResult < 0){
                continue;
            }
            else if(bufferedReadResult >= 0){ //소리를 읽은 경우

                for(int index = 0 ; index < buffer.length ; index++){
                    transformed[index] = (double)buffer[index];  //short -> double : findFrequency( )사용을 위해
                }

                double dom = findFrequency(transformed);
                //Log.d("dom : ", Double.toString(dom));


                if(startFlag && match(dom, HANDSHAKE_END_HZ)){
                    Log.d("packet : ", packet.toString());
                    Double[] packetArray = packet.toArray(new Double[packet.size()]);
                    List<Integer> data = this.extract_packet(packetArray);
                    Log.d("List data : ", data.toString());
                    String result = "";
                    for(int index = 0 ; index < data.size() ; index++){
                        result += Character.toString((char)((int)data.get(index)));
                    }
                    Log.d("RESULT : ", result);

                    packet.clear();
                    startFlag = false;
                }
                else if(startFlag){
                    packet.add(dom);
                }
                else if(match(dom, HANDSHAKE_START_HZ)){
                    startFlag = true;
                }

            }
        }

    }

    private boolean match(double freq1, int freq2){
        return abs(freq1 - freq2) < 20;
    }

    private List<Integer> extract_packet(Double[] packet) {
        int[] bitChunks = new int[packet.length / 2];
        int[] chunks = new int[packet.length / 2 - 1];
        double[] freqs = new double[packet.length / 2];

        for(int index = 0, paramIndex = 0 ; index < packet.length / 2 ; index++, paramIndex = paramIndex + 2){
            freqs[index] = packet[paramIndex];
        }

        for(int index = 0 ; index < freqs.length ; index++){
            bitChunks[index] = (int)Math.round((freqs[index] - START_HZ) / STEP_HZ);
        }

        for(int index = 1 ; index < freqs.length ; index++){
            if(bitChunks[index] >= 0 && bitChunks[index] < 16){
                chunks[index - 1] = bitChunks[index];
            }
        }

        List<Integer> result = decodeBitchunks(chunks);

        return result;
    }

    private List<Integer> decodeBitchunks(int[] chunks){
        List<Integer> outBytes = new ArrayList<>();

        int nextReadChunk = 0;
        int nextReadBit = 0;

        int value = 0; //byte
        int bitsLeft = 8;

        while(nextReadChunk < chunks.length){
            int canFill = 4 - nextReadBit;
            int toFill = min(bitsLeft, canFill);
            int offSet = 4 - nextReadBit - toFill;
            value = value << toFill;
            int shifted = chunks[nextReadChunk] & (((1 << toFill) - 1) << offSet);
            value = value | (shifted >> offSet);
            bitsLeft -= toFill;
            nextReadBit += toFill;

            if(bitsLeft <= 0){
                outBytes.add(value);
                value = 0;
                bitsLeft = 8;
            }

            if(nextReadBit >= 4){
                nextReadChunk += 1;
                nextReadBit -= 4;
            }
        }

        return outBytes;
    }
}

