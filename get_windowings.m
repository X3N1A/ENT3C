function [BOX_SIZES] = get_windowings(Resolution)

if Resolution==5e3
    BOX_SIZES_S=100;BOX_SIZES_E=250;
elseif Resolution==10e3
    BOX_SIZES_S=75;BOX_SIZES_E=200;
elseif Resolution==25e3
    BOX_SIZES_S=50;BOX_SIZES_E=150;
elseif Resolution==40e3
    BOX_SIZES_S=15;BOX_SIZES_E=200;    
elseif Resolution==50e3
    BOX_SIZES_S=20;BOX_SIZES_E=100;
elseif Resolution==100e3
    BOX_SIZES_S=15;BOX_SIZES_E=70;
elseif Resolution==250e3
    BOX_SIZES_S=10;BOX_SIZES_E=50;
elseif Resolution==500e3
    BOX_SIZES_S=10;BOX_SIZES_E=40;
elseif Resolution==1e6
    BOX_SIZES_S=10;BOX_SIZES_E=20;
end
    BOX_SIZES=round(linspace(BOX_SIZES_S,BOX_SIZES_E,9));
end

% for BOX_SIZE=BOX_SIZEs
%     WN=WNs(k);
%     WS=ceil((size(microC,1)-BOX_SIZE)/(WN-1));
%     WN=1+floor((size(microC,1)-BOX_SIZE)./WS);
%     [Resolution size(microC,1) BOX_SIZE  WN WS]
% end