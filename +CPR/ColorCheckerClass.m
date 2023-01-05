%% ColorChecker superclass in 2D and AR/VR
classdef ColorCheckerClass < handle
    %COLORCHECKERCLASS ColorChecker superclass
    %   See also ColorGauge, RezChecker, ColorCheckerClass_1x24

    properties
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

        target_type = 1;                      % 1: ColorChecker, 2: ColorGauge, 3: RezChecker
        patch_n = 24;                         % total number of patches including empty ones 
        patch_col = 6;                        % layout of the original target
        patch_row = 4;                        % layout of the original target
        all_patches = [1:24];                 % all patches that not empty
        grayscale_patches = [19:24];          % list of gray patches
        chromatic_patches = [1:18];           % list of chromatic patches
        white_patch = 19;                     % the whitest patch 
        
        dx = 120;                             % size used in the old order_1D code
    end

    methods

        function im = getImage_with_index (obj,colorchecker_rgb,idx)
            % Get an image of patches of idx

            if nargin <= 1
                colorchecker_rgb = obj.colorchecker;
                idx = 1:obj.patch_n;
            end

            wx = min(obj.patch_col,size(colorchecker_rgb,1));
            wy = obj.patch_row;
            dx = obj.dx;

            im = obj.getImage(colorchecker_rgb);

            for i=0:size(colorchecker_rgb,1)-1
                y=floor(i/wx);
                x=mod(i,wx);

                posx = x*dx + 10;
                posy = y*dx + 10;

                position = [posx posy];
                text = sprintf('%02d',idx(i+1));
                im = insertText(im,position,text,'FontSize',28);

            end
        end

        function im = getImage (obj,colorchecker_rgb)
            % Get an image of the input RGB values

            if nargin <= 1
                colorchecker_rgb = obj.colorchecker;
            end

            wx = min(obj.patch_col,size(colorchecker_rgb,1));
            wy = obj.patch_row;
            dx = obj.dx;

            im = zeros(dx*wy,dx*wx,3);

            for i=0:size(colorchecker_rgb,1)-1
                y=floor(i/wx);
                x=mod(i,wx);
                for channel=1:3
                    im(dx*y+1:dx*(y+1),dx*x+1:dx*(x+1),channel)=colorchecker_rgb(i+1,channel)/255;
                end;
            end

            %             imshow(im);
            %             axis equal;

            % imwrite(a,'color24.tif');
        end

        function arvr (obj,colorchecker_rgb)
            % Generate HTML for A-frame

            nargin
            if nargin <= 1
                colorchecker_rgb = obj.colorchecker;
            end

            colorchecker_lab = rgb2lab(colorchecker_rgb);

            fileID = fopen('colorchecker_in_vr.html','w');


            fprintf(fileID,'%s','<html>  <head>    <script src="https://aframe.io/releases/1.2.0/aframe.min.js"></script>  </head>  <body>    <a-scene>');

            for i = 1:length(colorchecker_lab)

                % coordinate for the color tree
                % 40 pages dividing 360 degrees
                lab = colorchecker_lab(i,:);

                sc = 0.00010;
                y = lab(1)*sc + 0.5;
                x = lab(2)*sc;
                z = lab(3)*sc;
                r = 0.030;

                % create html A-frame
                rgbuint8 = uint8(colorchecker_rgb(i,:));
                rgbcode = sprintf('%02X%02X%02X',rgbuint8(1),rgbuint8(2),rgbuint8(3));

                % output
                if nargin <= 1
                    fprintf(fileID,'<a-sphere position="%.4f %.4f %.4f" radius="%.4f" color="#%s"></a-sphere>\n',x,y,z,r,rgbcode);
                else
                    fprintf(fileID,'<a-box position="%.4f %.4f %.4f" depth="%.4f" height="%.4f" width="%.4f" color="#%s"></a-box>\n',x,y,z,r,r,r,rgbcode);
                end

            end

            fprintf(fileID,'%s','<a-sky color="#ECECEC"></a-sky>   </a-scene>  </body></html>            ');

            fclose(fileID);
        end

        function arvrx2 (obj,rgb_reference,rgb_device,html_fn)
            % Generate HTML for A-frame
            % usage: cc.arvrx2(cpr.rgbOriginal,cpr.rgbReproduced)
            % must input 2 datasets of the same size
            assert(nargin >= 3);

            % lab
            lab_reference = rgb2lab(rgb_reference);
            lab_device = rgb2lab(rgb_device);

            if nargin < 4
                html_fn = 'colorchecker_in_vr.html';
            end

            fileID = fopen(html_fn,'w');

            %
            % before
            %
            fprintf(fileID,'%s\n','<html>  <head>    <script src="https://aframe.io/releases/1.2.0/aframe.min.js"></script>  </head>  <body>    <a-scene>');

            %
            % reference
            %
            for i = 1:length(lab_reference)

                % coordinate
                xyz_reference = lab2xyz(lab_reference(i,:));
                xyz_device = lab2xyz(lab_device(i,:));
                
                % ball size
                r = 0.030;

                % create html A-frame
                rgbcode_reference = rgb2code(rgb_reference(i,:));
                rgbcode_device = rgb2code(rgb_device(i,:));

                %
                % output
                %

                % balls
                fprintf(fileID,'<a-sphere position="%.4f %.4f %.4f" radius="%.4f" color="#%s"></a-sphere>\n', ...
                    xyz_reference(1),xyz_reference(2),xyz_reference(3),r,rgbcode_reference);
                
                % boxes
                fprintf(fileID,'<a-box position="%.4f %.4f %.4f" depth="%.4f" height="%.4f" width="%.4f" color="#%s"></a-box>\n', ...
                    xyz_device(1),xyz_device(2),xyz_device(3),r,r,r,rgbcode_device);
                
                % lines
                fprintf(fileID,'<a-entity line="start: %.4f %.4f %.4f; end: %.4f %.4f %.4f; color: #%s"></a-entity>\n', ...
                    xyz_reference(1),xyz_reference(2),xyz_reference(3), ...
                    xyz_device(1),xyz_device(2),xyz_device(3), ...
                    rgbcode_reference);

            end

            %
            % after
            %
            fprintf(fileID,'%s\n','<a-sky color="#ECECEC"></a-sky>   </a-scene>  </body></html>            ');

            fclose(fileID);

            return

            function xyz = lab2xyz (lab)
                sc = 0.00010;
                y = lab(1)*sc + 0.5;
                x = lab(2)*sc;
                z = lab(3)*sc;
                xyz = [x y z];
            end
            
            function code = rgb2code (rgb)
                rgbuint8 = uint8(rgb);
                code = sprintf('%02X%02X%02X',rgbuint8(1),rgbuint8(2),rgbuint8(3));
            end
        end

    end
end

