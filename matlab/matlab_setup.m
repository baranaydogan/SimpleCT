% ------------------------------------------------------------
% Dogu Baran Aydogan - baran.aydogan@gmail.com
% 21.07.2014
% ------------------------------------------------------------

% Compiles SimpleCT's Matlab interface
% Only tested in Linux


archstr = computer('arch');
disp(' ');
if archstr(1) == 'g'
    % Linux - tested
    disp(['Compiling SimpleCT Matlab interface for Linux (', archstr, ')...']);
    copyfile('../src/SimpleCT.cpp', '.');
    copyfile('../src/SimpleCT.h', '.');
    mex ContourTree.cpp SimpleCT.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
    delete('SimpleCT.cpp');
    delete('SimpleCT.h');
    disp('Done');
elseif archstr(1) == 'w'
    % Windows - NOT TESTED
    disp(['Compiling SimpleCT Matlab interface for Windows (', archstr, ')']);
    copyfile('..\src\SimpleCT.cpp', '.');
    copyfile('..\src\SimpleCT.h', '.');
    mex ContourTree.cpp SimpleCT.cpp COMPFLAGS="$COMPFLAGS /openmp"
    delete('SimpleCT.cpp');
    delete('SimpleCT.h');
    disp('Done');
elseif archstr(1) == 'm'
    % Mac - NOT TESTED
    disp(['Compiling SimpleCT Matlab interface for Mac (', archstr, ')']);
    copyfile('../src/SimpleCT.cpp', '.');
    copyfile('../src/SimpleCT.h', '.');
    mex ContourTree.cpp SimpleCT.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
    delete('SimpleCT.cpp');
    delete('SimpleCT.h');
    disp('Done');
end

example