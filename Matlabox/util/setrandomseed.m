function [stream, prevstream] = setrandomseed(rseed)

    prevstream = RandStream.getGlobalStream();
    
    stream = RandStream('mt19937ar', 'Seed', rseed);
    RandStream.setGlobalStream(stream);
    