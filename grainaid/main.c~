/* Julius Raskevicius, GrainAid, 2016-2021
    Description: polyphonic granular synthesis engine written with PortAudio library. 
    Takes a sound file as an input and produces a granulated audio output for 120 layers (see static parameters below).
    Works well as a granular stretcher (see Examples for audio output)
    Parameters: 
        filename - string, path to the input file
        start_sample - offset for sample start in seconds
        amp_factor - try adjusting this, if the sound is too quiet or loud (approx. 1/2 kLayerNumber)
*/

/*
 For further tutorials on building PortAudio plugins search for pa_skeleton.c file.
 
 It is a framework for building a C program that uses PortAudio
 loosely based on Phil Burk's paex_saw.c example
*/


/*
 * This program uses the PortAudio Portable Audio Library.
 * For more information see: http://www.portaudio.com
 * Copyright (c) 1999-2000 Ross Bencina and Phil Burk
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files
 * (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge,
 * publish, distribute, sublicense, and/or sell copies of the Software,
 * and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
 * ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * The text above constitutes the entire PortAudio license; however,
 * the PortAudio community also makes the following non-binding requests:
 *
 * Any person wishing to distribute modifications to the Software is
 * requested to send the modifications to the original developer so that
 * they can be incorporated into the canonical version. It is also
 * requested that these non-binding requests be included along with the
 * license above.
 */



#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <time.h>
#include "portaudio.h"
#include "sndfile.h"

// Controls:
#define kSampleRate   (44100)
#define kFramesPerBuffer (256)
#define kChan 2

#define kLayerNumber 100        // Number of layers
                                // CAUTION: some things are normalized by layer count, setting to low value without reducing kAmpFactor will produce loud audio!
//#define kAmpFactor 60           // Amp factor.
                                // CAUTION: large factor = ear damage / crash.

#define grainLength 44100 / 15  // Grain duration (samples). Grains are windowed, thus defines windowing effect too
#define kLocPerc 1              // Random deviation (%) of grain start
#define rRate 0.1               // Read rate (< 1 for slowdown, >>1 for randomness). Number of samples to progress for each sample in buffer
#define kGChance 1              // Chance of grain unmuted (~density, ~channel_separation).
//#define kStartSamp 44100 * 0    // Offset in buffer
#define offset 0.02             // A-S-R offset (A-S & S-R)

/* define a struct for  passing data to your callback function */
typedef struct {
	float *buf[2];                          // array holding pointers to L and R channels.
    double loc;                             // location of global seek head

    // "0 phase" grains
    double gLoc[kLayerNumber][kChan];       // location of specific layer/channel grain in buffer
    double gLen[kLayerNumber][kChan];       // grain length of specific layer/channel
    double gLeft[kLayerNumber][kChan];      // remaining samples to play in grain of specific layer/channel
    
    // "1/3 phase" grains
    double gPLoc[kLayerNumber][kChan];
    double gPLen[kLayerNumber][kChan];
    double gPLeft[kLayerNumber][kChan];

    // "2/3 phase" grains
    double gPPLoc[kLayerNumber][kChan];
    double gPPLen[kLayerNumber][kChan];
    double gPPLeft[kLayerNumber][kChan];
 
    int unMute[kLayerNumber][kChan];        // unmute status of layer/channel (1 = grain is silent)
    int unMuteP[kLayerNumber][kChan];
    int unMutePP[kLayerNumber][kChan];
    
    int amp;                                // amplification factor
   
    double items;                           // total number of samples across both channels of a sample
} myData;

double gLDev(double perc);

int makeGrain(double chance);

float *buf1, *buf2;

// sndfilelib declarations:
// ************************


SNDFILE*  sf_open    (const char *path, int mode, SF_INFO *sfinfo) ;
/* When opening a file for read, the format field should be set to zero before calling sf_open(). 
  The only exception to this is the case of RAW files where the caller has to set the samplerate, 
  channels and format fields to valid values. All other fields of the structure are filled in by the library.
*/

/*
 SFM_READ    - read only mode
 SFM_WRITE   - write only mode
 SFM_RDWR    - read/write mode
 */

/* SF_INFO struct as defined in sndfile.h
 
 typedef struct
 {    sf_count_t  frames ;
 int         samplerate ;
 int         channels ;
 int         format ;
 int         sections ;
 int         seekable ;
 } SF_INFO ;
 */

sf_count_t  sf_read_double  (SNDFILE *sndfile, double *ptr, sf_count_t items) ;

int ct = 0, l, ch;

double *buf;
double f, items, num, num2;
int sr, c;


/*******************************************************************/

/*  the CALLBACK ROUTINE FOR PRODUCING AUDIO */

/* This routine will be called by the PortAudio engine when audio is needed.
 ** It may called at interrupt level on some machines so don't do anything
 ** that could mess up the system like calling malloc() or free().
 */

static int PaTestCallback(const void *inputBuffer, void *outputBuffer,
                          unsigned long framesPerBuffer,
                          const PaStreamCallbackTimeInfo* timeInfo,
                          PaStreamCallbackFlags statusFlags,
                          void *userData ){
    
	myData *data = (myData*)userData;       // user-defined input data pointer
	float *out = (float*)outputBuffer;      // output buffer
    float samp[kChan] = {0, 0};             // read head for each channel
    int i, d, ip, dp;
    double env, envp, envpp;
    
    while(framesPerBuffer--) {
        samp[0] = 0;
        samp[1] = 0;

        for (l = 0; l < kLayerNumber; l++) {
            for (ch = 0; ch<2; ch++) {

                // 1. If grain playback is over, generate new grains
                // "1/3 and 2/3 phase" grains double the material with 1/3 grain phase shift to remove stutter
                if (data->gLeft[l][ch] == 0) {
                    
                    // New "0 phase" grain
                    if (makeGrain(kGChance))
                        data->unMute[l][ch] = 1;
                    else
                        data->unMute[l][ch] = 0;
                    
                    data->gLoc[l][ch] = (int) (data->loc * gLDev(kLocPerc)) % (int) data->items;
                    data->gLen[l][ch] = (int) grainLength; // * gLDev(5));
                    data->gLeft[l][ch] = data->gLen[l][ch];

                }
                
                // New "1/3 phase" grain
                if (data->gPLeft[l][ch] == 0) {
                    if (makeGrain(kGChance)) {
                        data->unMuteP[l][ch] = 1;
                    }
                    else {
                        data->unMuteP[l][ch] = 0; 
                    }

                    data->gPLoc[l][ch] = (int) (data->gLoc[l][ch] * gLDev(kLocPerc)) % (int) data->items;
                    data->gPLen[l][ch] = data->gLen[l][ch];
                    data->gPLeft[l][ch] = data->gPLen[l][ch];
                }

                // New "2/3 phase" grain
                if (data->gPPLeft[l][ch] == 0) {
                    if (makeGrain(kGChance)) {
                        data->unMutePP[l][ch] = 1;
                    }
                    else {
                        data->unMutePP[l][ch] = 0; 
                    }

                    data->gPPLoc[l][ch] = (int) (data->gLoc[l][ch] * gLDev(kLocPerc)) % (int) data->items;
                    data->gPPLen[l][ch] = data->gLen[l][ch];
                    data->gPPLeft[l][ch] = data->gPPLen[l][ch];
                }
        
                // 2. Calculate grain envelopes (linear attack-sustain-decay)              
                // Env. "0 phase"
                double r = (1 - data->gLeft[l][ch] / data->gLen[l][ch]);  // ramp
                if (r < offset) {
                        env = r / offset;
                }
		        else if (r >= (1 - offset)) {
			        env = (1 - r) / offset;;
                }
		        else {
                    env = 1.;                
                }

                // Env. "1/3 phase"
                double r_p = (1 - data->gPLeft[l][ch] / data->gPLen[l][ch]);  // phase ramp
                if (r_p < offset) {
                    envp = r_p / offset;
                }

		        else if (r_p >= (1 - offset)) {
			        envp = (1 - r_p) / offset;
                }
		        else {
                    envp = 1.;                 
                }

                // Env. "2/3 phase"
                double r_pp = (1 - data->gPPLeft[l][ch] / data->gPPLen[l][ch]);  // phase ramp
                if (r_pp < offset) {
                    envpp = r_pp / offset;
                }

		        else if (r_pp >= (1 - offset)) {
			        envpp = (1 - r_pp) / offset;
                }
		        else {
                    envpp = 1.;                 
                }

                // TODO: DEBUG envelopes
//                if ((l == 0) && (ch == 0)) {
//                    printf("%f, %f, %f, %f, %f, %f \n", env, envp, envpp, r, r_p, r_pp);
//                }
//          
                // 3. Calculate output
                samp[ch] += data->buf[ch][(int) data->gLoc[l][ch]] * env /2.0 * data->unMute[l][ch];
                if (data->gPLoc[l][ch] < data->items) {
                    samp[ch] += data->buf[ch][(int) data->gPLoc[l][ch]] * envp /2.0 * data->unMuteP[l][ch];
                }
                if (data->gPPLoc[l][ch] < data->items) {
                    samp[ch] += data->buf[ch][(int) data->gPPLoc[l][ch]] * envpp /2.0 * data->unMutePP[l][ch];
                }
//                samp[ch] /= kLayerNumber/kAmpFactor;  // reduce overall level
                samp[ch] /= kLayerNumber/data->amp;  // reduce overall level

                // 4. Update grain position counters;
                data->gLoc[l][ch] =  data->gLoc[l][ch] + 1;
                data->gLeft[l][ch] = data->gLeft[l][ch] - 1;
                
                data->gPLoc[l][ch] =  data->gPLoc[l][ch] + 1;
                data->gPLeft[l][ch] = data->gPLeft[l][ch] - 1;

                data->gPPLoc[l][ch] =  data->gPPLoc[l][ch] + 1;
                data->gPPLeft[l][ch] = data->gPPLeft[l][ch] - 1;
            }
        }
                
        // 5. Once we add each sample from each grain, progress the global seek head by rRate
        data->loc += rRate;
        if (data->loc > data->items) {
            data->loc = 0;
        }
        
        // 6. Write calculated sample to the output buffer
        for (ch = 0; ch < 2; ch++) {
            *out++ = samp[ch];
        }
    }
	return 0;
}

/*******************************************************************/

int main(int argc, char **argv)
{
	PaStream *stream;  /* declare a stream variable for PortAudio */
	PaError err;  /* declare an error variable */
	myData data; // declare structure for user data (if needed)

    float *stereoBuf;
    int ctr = 0;
    int ctrL = 0;
    int ctrR = 0;
    
    srand((int) time(NULL));

    // Get user input from command line
    char* filePath = argv[1];           // filename
    char* startSample_str = argv[2];    // starting point (location in sample (s))
    char* amp_str = argv[3];            // amplification factor
    char *null_ptr;

    long startSample = 0;
    long amp = 0;

    startSample = strtol(startSample_str, &null_ptr, 10);
    startSample = startSample * 44100;  // convert input to seconds
    amp = strtol(amp_str, &null_ptr, 10);
    printf("Starting point (s): %ld\n", startSample);
    printf("Amplification factor: %ld\n\n", amp);

    // Open file with libsndfile
    /////////////////////////////
    
    SF_INFO mSfInfo;
    mSfInfo.format = 0;
    
//    SNDFILE *sf = sf_open(kFN, SFM_READ, &mSfInfo);
    SNDFILE *sf = sf_open(filePath, SFM_READ, &mSfInfo);
    if (sf == NULL) {
        printf("Could't open file! \n");
        exit(0);
    }
    
    f = mSfInfo.frames;
    sr = mSfInfo.samplerate;
    c = mSfInfo.channels;
    items = c * f;
    data.items = items; // Items is frames * channels!
    
    // Read file into buffer
    stereoBuf = (float *) malloc(items*sizeof(float));
    num = sf_read_float(sf,stereoBuf,items);
    sf_close(sf);

    // Separate into L/R arrays
    data.buf[0] = (float *) malloc(f*sizeof(float));
    data.buf[1] = (float *) malloc(f*sizeof(float));
    
    for (ctr = 0; ctr<f; ctr++) {
        if (ctr % 2 == 0) {
            data.buf[0][ctrL++] = stereoBuf[ctr];
        }
        else {
            data.buf[1][ctrR++] = stereoBuf[ctr];
        }
    }

    /* Initialize data for use by callback. */
//    data.loc = kStartSamp;
    data.loc = startSample; // start location for main position is first sample buffer.
    
    for (l = 0; l < kLayerNumber; l++) {
        for (ch = 0; ch < 2; ch++) {
            data.gLoc[l][ch] = data.loc * gLDev(kLocPerc); // gLDev return range [1; 2]
            data.gLen[l][ch] = (int) (grainLength); // * gLDev(5)); // Cast to int, otherwise never 0.
            data.gLeft[l][ch] = 0;
            data.gPLeft[l][ch] = data.gLen[l][ch] * 1.0 / 3.0; // We offset the first "1/3 phase" grain
            data.gPLen[l][ch] = (int) (grainLength); 
            
            data.gPPLeft[l][ch] = data.gLen[l][ch] * 2.0 / 3.0; // We offset the first "2/3 phase" grain
            data.gPPLen[l][ch] = (int) (grainLength); 
            data.unMute[l][ch] = 0;
            data.unMuteP[l][ch] = 0;
            data.unMutePP[l][ch] = 0;
            data.amp = amp;
        }
    }

    // Display file data to user
    printf("File info: \n");
    printf("frames=%d\n",(int) f);
    printf("samplerate=%d\n",sr);
    printf("channels=%d\n", c);
    printf("items=%d\n\n", (int) items);
    
    // PortAudio...
    printf("PortAudio Test\n");  /* announce what is happening */
    
    /* Initialize library before making any other calls. */
	err = Pa_Initialize();
	if( err != paNoError ) {
        goto error;
    }
    
    /* Open an audio I/O stream. */
	err = Pa_OpenDefaultStream(
                               &stream,
                               0,              /* no input channels */
                               2,              /* stereo output */
                               paFloat32,      /* 32 bit floating point output */
                               kSampleRate,
                               kFramesPerBuffer,            /* frames per buffer */
                               PaTestCallback,  /* pointer to the callback function */
                               &data);          /* data we pass to the callback  */
	if( err != paNoError ) {
        goto error;
    }
    
	err = Pa_StartStream( stream );  /* start stream */
	if( err != paNoError ) {
        goto error;
    }
    
    /* Sit and wait for the user to type a Return. This gives PortAudio time to
     call the callback function and produce some audio for a while*/
    printf("Hit RETURN to stop.\n");
    getchar();
    
	err = Pa_StopStream( stream );  /* after RETURN has been hit, stop stream */
	if( err != paNoError ) {
        goto error;
    }
    
	err = Pa_CloseStream( stream );  /* close stream */
	if( err != paNoError ) {
        goto error;
    }
    
	Pa_Terminate(); /* terminate PortAudio */
	printf("Test finished.\n");
	return err;  /* normal return for this function */
    
error:
	Pa_Terminate();    /* If an error occured, terminate and then print the error message */
	fprintf( stderr, "An error occured while using the portaudio stream\n" );
	fprintf( stderr, "Error number: %d\n", err );
	fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
	return err;
}


double gLDev(double perc) {
    /* Calculate grain location deviation (grain jitter)
    @params:
        int perc (percentage shift)
    @return:
        double in range [0.75; 1.75], multiplier for grain location    
    */
     return 0.75 + (((double) perc / 100) * ((double) rand() / RAND_MAX));
}

int makeGrain(double chance) {
    /* Decide if a grain should be made given chance
    @params:
        double chance (prob. of grain audible)
    @return:
        int 0/1, boolean that determines if the grain is audible
    */
    float val;

    val =  ((double) rand() / RAND_MAX);
    return val < chance;
}


