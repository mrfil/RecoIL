function [X, dXdx, dXdy] = bump(N,w,xloc,yloc,opt,fmr)
   %
   %function X = bump(N,w,xloc,yloc,opt,fmr)
   %  this function is used for making blobs in images
   %    can be used to make field maps, etc.
   %     opt == 1, then do gaussian
   %     opt == 2, then do fermi
   %   xloc, yloc -> location of center
   %   w is width of the gaussian bump
   %        width of the flat region of fermi
   %   fmr is the transition quickness of the fermi filter,
   %         gives transition region width
   %Brad Sutton, U of Michigan, Sep. 2002
   
   [xvals,yvals] = meshgrid(1:N);
   
   
   if (opt == 1)
      d = ((xvals-xloc).^2+(yvals-yloc).^2);
      X = exp(-1./2.*(d./(w^2)));
      if nargout > 1
         dXdx = - (xvals - xloc).*X./(w.^2);
         dXdy = - (yvals - yloc).*X./(w.^2);
      end
   end
   
   
   if (opt ==2)
      d = sqrt((xvals-xloc).^2+(yvals-yloc).^2);
      X = 1./(1+exp((d-w)./fmr));
      if nargout > 1
         dXdx = - (xvals - xloc).*exp((d-w)./fmr)./(fmr.*d.*(exp((d-w)./fmr)+1).^2);
         dXdy = - (yvals - yloc).*exp((d-w)./fmr)./(fmr.*d.*(exp((d-w)./fmr)+1).^2);
      end
   end
