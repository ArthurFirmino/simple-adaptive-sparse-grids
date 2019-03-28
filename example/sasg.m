function obj = sasg(ndim,nlvl)
    % Create subfolder, add to path
    if exist('sasg_builds','dir') ~= 7 
        mkdir('sasg_builds')
    end
    addpath('./C___class_interface')
    addpath('./sasg_builds')
    % Check existence of desired mex, else compile
    if exist('sasg_builds/sasg_mex_d'+string(ndim), 'file') ~= 3
        mex('-IC:./C___class_interface','-IC:../',...
            'COMPFLAGS=$COMPFLAGS -O2 -DNDIMENSIONS='...
            +string(ndim),'-outdir','sasg_builds',...
            '-output','sasg_mex_d'+string(ndim),'sasg_mex.cpp');
    end
    obj = mex_interface(str2fun(char('sasg_mex_d'+string(ndim))), nlvl);
end