function img = mergePowerGridFileOutput(NSlices,NReps,NAvgs,NEchoes,NPhases,directory)
% mergePowerGridFileOutput(NSlices,NReps,NAvgs,NEchoes,NPhases)
% Forms complex N-dimensional images from the output of PowerGrid. Does not
% handle multiband reshuffling if that is required.

if nargin < 6
    directory = pwd;
end

for phs = 1:NPhases
    for eco = 1:NEchoes
        for avg = 1:NAvgs
            for rep = 1:NReps
                for slc = 1:NSlices
                    tmp1 = load_untouch_nii(sprintf([directory '/' 'img_Slice%i_Rep%i_Avg%i_Echo%i_Phase%i_mag.nii'],slc-1,rep-1,avg-1,eco-1,phs-1),1,1,1,1);
                    tmp2 = load_untouch_nii(sprintf([directory '/' 'img_Slice%i_Rep%i_Avg%i_Echo%i_Phase%i_phs.nii'],slc-1,rep-1,avg-1,eco-1,phs-1),1,1,1,1);
                    img(:,:,:,slc,rep,avg,eco,phs) = tmp1.img.*exp(1i*tmp2.img);
                end
            end
        end
    end
end

