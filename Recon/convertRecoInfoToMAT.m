function convertRecoInfoToMAT(dirname,rInfo,sen, FM, PhaseMaps)

if exist(dirname,'dir')
    error(['Directory ' dirname ' already exists.  Please remove first'])
    %delete(filename);
else
    mkdir(dirname);
    cd(dirname);
end

sen = col(sen);
save('sen.mat',col(sen));


FM = col(FM);
save('FM.mat',col(FM));


% if exists append PMaps to file
if nargin > 4
    PhaseMaps = col(PhaseMaps);
    save('PMaps.mat',PMaps);
end


% ASSUME WE HAVE ONLY ONE IMAGE FOR NOW

eco = 1;
phs = 1;
avg = 1;

for rep = 1:1
    for slc = 1:1
        for par = 1:rInfo.nPartitions
            for shot = 1:rInfo.nShots
                              
                % fill the data
                dataTemp = squeeze(rInfo.dataRead(shot,par,slc,avg,phs,eco,rep));
                % attach the trajectory
                kxTemp = col(rInfo.kx(rInfo.dataMask,shot,par,slc,avg,phs,eco,rep));
                kyTemp = col(rInfo.ky(rInfo.dataMask,shot,par,slc,avg,phs,eco,rep));
                kzTemp = col(rInfo.kz(rInfo.dataMask,shot,par,slc,avg,phs,eco,rep));
                ttTemp = col(rInfo.timingVec(rInfo.dataMask,shot,par,slc,avg,phs,eco,rep));
                
               kx(:,par,shot) = kxTemp;
               ky(:,par,shot) = kyTemp;
               kz(:,par,shot) = kzTemp;
               tt(:,par,shot) = ttTemp;
               data(:,par,shot,:) = dataTemp;
               
            end % shot loop
        end % partition loop
    end % slice loop
end % rep loop

kx = col(kx);
ky = col(ky);
kz = col(kz);
tt = col(tt);
data = col(data);

save('kx.mat',kx);
save('ky.mat',ky);
save('kz.mat',kz);
save('tt.mat',tt);
save('data.mat',data);

end

