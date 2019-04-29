function mapout = resample_map_resolution(mapin,Nout,Noutsl,fov_in, fov_out)
   % function mapout = resample_map_resolution(mapin,recoInfo.N,recoInfo.nsl,fov_in, fov_out)
   % This function resamples a map to match a target resolution
   % WARNING: IT ASSUMES THAT THE FOV in slice direction is the same between input
   %     output.
   % mapin may be complex
   %
   % Last edited: 8/16/2013. BPS: added fov_in, fov_out
   
   
   if isreal(mapin)
      flag_phase = 0;
   else
      flag_phase = 1;
   end
   
   nsl = size(mapin,3);
   nszxy = size(mapin);
   nx = nszxy(2);
   
   
   if ~exist('Noutsl','var')
      Noutsl = nsl;
   end
   if ~exist('fov_in','var')
      fov_in = 240;
   end
   if ~exist('fov_out','var')
      fov_out = 240;
   end
   
   
   
   
   if ~(nszxy(1) == nszxy(2))
      sprintf('resampling only works for square inputs')
   end
   
   input_vec = (([0:nx-1] - nx/2 + 1/2)./nx)*fov_in;
   output_vec = (([0:Nout-1] - Nout/2 + 1/2)./Nout)*fov_out;
   
   %mapout_mag = zeros(Nout,Nout,nsl);
   %mapout_phs = zeros(Nout,Nout,nsl);
   
   for sliceIndex = 1:nsl
      if ~flag_phase
         mapout_mag(:,:,sliceIndex) = interp2(input_vec,input_vec',real(mapin(:,:,sliceIndex)),output_vec,output_vec','cubic');
      else
         mapout_mag(:,:,sliceIndex) = interp2(input_vec,input_vec',abs(mapin(:,:,sliceIndex)),output_vec,output_vec','cubic');
         mapout_phs(:,:,sliceIndex) = interp2(input_vec,input_vec',angle(mapin(:,:,sliceIndex)),output_vec,output_vec','nearest');
      end
   end
   
   if flag_phase
      mapout = mapout_mag.*exp(1i*mapout_phs);
   else
      mapout = mapout_mag;
   end
   
   mapout(isnan(mapout)) = 0;
   
   
   
   if ~(nsl == Noutsl)
      mapout2_mag = zeros(size(mapout,1),size(mapout,2),Noutsl);
      mapout2_phs = zeros(size(mapout,1),size(mapout,2),Noutsl);
      
      nslin_vec = ([-nsl/2:nsl/2-1]+0.5)./nsl;
      nslout_vec = ([-Noutsl/2:Noutsl/2-1]+0.5)./Noutsl;
      
      for ii = 1:size(mapout,1)
         for jj = 1:size(mapout,2)
            if flag_phase
               mapout2_phs(ii,jj,:) = interp1(nslin_vec(:),squeeze(unwrap(angle(mapout(ii,jj,:)))),nslout_vec(:));
               mapout2_mag(ii,jj,:) = interp1(nslin_vec(:),squeeze(abs(mapout(ii,jj,:))),nslout_vec(:));
            else
               mapout2_mag(ii,jj,:) = interp1(nslin_vec(:),squeeze(real(mapout(ii,jj,:))),nslout_vec(:));
            end
            
         end
      end
      if flag_phase
         mapout2 = mapout2_mag.*exp(1i*mapout2_phs);
      else
         mapout2 = mapout2_mag;
      end
      
      % FIX INTERPOLATION AT EDGES
      if (nslin_vec(1) >= nslout_vec(1))
         slvec_comp = nslout_vec - nslin_vec(1);
         beg_sl = find(slvec_comp>0);
         sprintf('Fixing beginning slice due to no overlap at front edge \n')
         for ii = 1:(beg_sl(1)-1)
            mapout2(:,:,ii) = mapout2(:,:,beg_sl(1));
         end
      end
      if (nslin_vec(end) <= nslout_vec(end))
         sprintf('Fixing end slice due to no overlap at back edge \n')
         slvec_comp = nslout_vec - nslin_vec(end);
         end_sl = find(slvec_comp<0);
         for ii = (end_sl(end)+1):Noutsl
            mapout2(:,:,ii) = mapout2(:,:,end_sl(end));
         end
      end
      
      mapout = mapout2;
   end
   
   mapout(isnan(mapout)) = 0;
   
   
   
