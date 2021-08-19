# GrainAid
    * Polyphonic granular synthesis engine. Takes a sound file as an input and produces a granulated audio output(see static parameters below)
    * Works well as a granular stretcher (see /examples for audio output)
    * Written in C with PortAudio and libsndfile

## Building:
    Linux:
    1) Install dependencies:
       # libsndfile - library for reading and writing sound files
       sudo apt-get install lsndfile1

       # PortAudio - library for audio playback and recording (needed to define the playback callback)
       sudo apt-get install libasound-dev portaudio19-dev libportaudio2 libportaudiocpp0
       sudo apt-get install ffmpeg libav-tools

    2) clang main.c -lportaudio -lsndfile -o granulator.out
    MacOS:
        See the XCodeProject in the repository

## Usage
Start from the command line. Compiled binary is located in /bin (for Linux). 

ATTENTION: There are no checks, if you input a parameter that is not valid, expect harsh noise.
Parameters:
        filename - string, path to the input file (has to be stereo file)
        start_sample - offset for sample start in seconds
        amp_factor - try adjusting this, if the sound is too quiet or loud (approx. 1/2 kLayerNumber, so in the range of 50-60)

Examples:
```
$ ./granulator ../samples/chimes.wav 1 50
$ ./granulator ../samples/orchestra.wav 0 50
$ ./granulator ../samples/dog.wav 0 50
```

## Issues:
    "...underrun occurred" - too much load for the CPU, reduce kLayerNumber and recompile 
    <HARSH NOISE> - check if the input parameters are of the right type and in the right range
