%% Superclass for all color test targets
% 12/1/2023

% called (including subclasses ColorGauge, RezChecker, and SFRPlus) by cprapp with ct:
%   getImage()
%   patch_n

classdef ColorCheckerClass < handle
    %COLORCHECKERCLASS ColorChecker superclass
    %   Default is the Gretag Macbeth ColorChecker
    %   See also ColorGauge, RezChecker, ColorCheckerClass_1x24

    %%
    properties
        %%

        name = 'ColorChecker'              % name of the test target

        % default RGB values of each patch
        colorchecker = [
            115 82  68;
            194 150 130;
            98  122 157;
            87  108 67;
            133 128 177;
            103 189 170;

            214 126 44;
            80  91  166;
            193 90  99;
            94  60  108;
            157 188 64;
            224 163 46;

            56  61  150;
            70  148 73;
            175 54  60;
            231 199 31;
            187 86  149;
            8   133 161;

            243 243 242;
            200 200 200;
            160 160 160;
            122 122 121;
            85  85  85;
            52  52  52
            ];

        target_type = 1;                      % 1: ColorChecker, 2: ColorGauge, 3: RezChecker, 4: SFRPlus
        patch_n = 24;                         % total number of patches including the empty ones

        patch_col = 6;                        % layout of the original target - number of columns
        patch_row = 4;                        % layout of the original target - number of rows
        all_patches = [1:24];                 % all patches that are not empty
        grayscale_patches = [19:24];          % list of the grayscale patches
        chromatic_patches = [1:18];           % list of the chromatic patches
        white_patch = 19;                     % the whitest patch

        dx = 120;                             % size used in the old order_1D code

    end

    %%
    methods

        %%
        function im = getImage_with_index (obj, colorchecker_rgb, idx)
            % Get an annotated image of the input RGB values
            % the template is filled with the input values

            % if colorchecker_rgb is absent, use the default values
            if nargin <= 1
                colorchecker_rgb = obj.colorchecker;
                idx = 1:obj.patch_n;
            end

            wx = min(obj.patch_col,size(colorchecker_rgb,1));          % how many columns
            wy = obj.patch_row;                                        % how many rows
            dimx = obj.dx;                                             % size of each patch in #pixel

            im = obj.getImage(colorchecker_rgb);

            for i=0:size(colorchecker_rgb,1)-1
                y = floor(i/wx);
                x = mod(i,wx);

                posx = x*dimx + 10;
                posy = y*dimx + 10;

                % need computer vision toolbox!
                position = [posx posy];
                text = sprintf('%02d',idx(i+1));
                im = insertText(im,position,text,'FontSize',28);
            end
        end

        %%
        function im = getImage (obj, colorchecker_rgb)
            % Get an image of the input RGB values
            % the template is filled with the input values

            % if colorchecker_rgb is absent, use the default values
            if nargin <= 1
                colorchecker_rgb = obj.colorchecker;
            end

            wx = min(obj.patch_col,size(colorchecker_rgb,1));         % how many columns
            wy = obj.patch_row;                                       % how many rows
            dimx = obj.dx;                                            % size of each patch in #pixel  

            im = zeros(dimx*wy,dimx*wx,3);                            % construct the image

            for i = 0:size(colorchecker_rgb,1)-1                        % for each patch 
                y = floor(i/wx);                                           % y location
                x = mod(i,wx);                                             % x location
                for channel=1:3                                             % assign R/G/B   
                    im(dimx*y+1:dimx*(y+1),dimx*x+1:dimx*(x+1),channel) = double(colorchecker_rgb(i+1,channel))/255;
                end
            end
        end

    end

end

