%% Colorimetrical data (CIEXYZ or CIELAB) of ColorChecker
classdef DeviceData < handle
    %DEVICEDATA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % source: 'spd', 'XYZ', or 'LAB'
        source = ''
        
        Lab = zeros(24,3)
        spec
        XYZ
        XYZ_n
        rgb
    end
    
    methods
        
        function obj = DeviceData
        end
        
        function construct_by_CL500a (obj,cl500a_xls_filepath)
            %%CONSTRUCT_BY_CL500a Import spectral data from CL500A via Excel

            % C:\Users\Wuchihlei\Documents\GitHub\MIS_acquire\@LightSim\cl500a_0608\CL500a_nec.xlsx
            spec_360_780 = xlsread(cl500a_xls_filepath);
            spec_380_780 = spec_360_780(:,21:end);
            
            obj.source = 'spd';
            obj.spec = spec_380_780;
            
            obj.spec2lab;
            
        end
        
        function construct_by_cs2000 (obj,cs2000_mat_filepath)
            %%CONSTRUCT_BY_CS2000 Import spectral data from CS2000 via .mat
            % during light matching

            % C:\Users\Wuchihlei\Documents\GitHub\MIS_acquire\@LightSim\cl500a_0608\colorchecker_test_result_nec.mat
            load(cs2000_mat_filepath,'spd_ol490_24')
            
            spec = zeros(24,401);
            for i = 1:24
                spec(i,:) = spd_ol490_24{i};
            end
            
            obj.source = 'spd';
            obj.spec = spec;
            
            obj.spec2lab;
        end
        
        function construct_by_target (obj,cs2000_mat_filepath)
            %%CONSTRUCT_BY_TARGET Import spectral data from target via .mat
            % during light matching

            % C:\Users\Wuchihlei\Documents\GitHub\MIS_acquire\@LightSim\cl500a_0608\colorchecker_test_result_nec.mat
            load(cs2000_mat_filepath,'spd_disp_24')
            
            spec = zeros(24,401);
            for i = 1:24
                spec(i,:) = spd_disp_24{i};
            end
            
            obj.source = 'spd';
            obj.spec = spec;
            
            obj.spec2lab;
        end
        
        function construct_by_i1 (obj,i1_mat_filepath)
            %%CONSTRUCT_BY_i1 Import spectral data from i1 via Excel
            % exported from i1Match and extracted to .mat

            % C:\Users\Wuchihlei\Documents\GitHub\MIS_acquire\@LightSim\nec_i1_0608\i1data_nec.mat
            load(i1_mat_filepath,'cxfdata')
            
            % cxfdata is 24x36 (patch x wavelength)
            % manually copy-n-paste from Excel to Matlab
            
            patch_order_in_cxf = [1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 23 24 3 4 5 6 7 8 9];
            spec_order = 380:10:730;
            
            patch_row_to_id = [[1:24]' patch_order_in_cxf'];
            patch_row_to_id_sorted = sortrows(patch_row_to_id,2);
            patch_id_to_row = [patch_row_to_id_sorted(:,2) patch_row_to_id_sorted(:,1)];
            
            spec_380_730 = cxfdata(patch_id_to_row(:,2),:);
            
            %
            % interpolate
            %
            spec_380_780 = zeros(24,401);
            for i = 1:24
                spec_i = interp1(spec_order,spec_380_730(i,:),380:780,'spline','extrap');
                spec_380_780(i,:) = spec_i;
            end
            
            obj.source = 'spd';
            obj.spec = spec_380_780;
            
            obj.spec2lab;
        end
        
        function construct_by_spyder (obj,i1_mat_filepath, mode)
            %%CONSTRUCT_BY_TARGET Import spectral data from .mat
            % manually entered in Excel and extracted to .mat
            % choose one of the 4 backlight modes
            
            % mode = 3
            modename = {'White LED','Standard LED','General','GB LED'};
            
            load(i1_mat_filepath,'spyder')
            
            cc = ColorConversionClass;
            
            colstart = (mode-1)*3 + 1;
            
            XYZ_cck = spyder(:,colstart:colstart+2);
            LAB_cck = CPR.ColorClass.XYZ2lab(XYZ_cck,XYZ_cck(19,:));
            
            obj.source = 'XYZ';
            obj.Lab = LAB_cck;
            
        end
        
        function spec2lab (obj)
            %%SPEC2LAB Convert spectrum to CIELAB

            for i = 1:24
                spec_i = obj.spec(i,1:10:end)';
                spec_white = obj.spec(19,1:10:end)';
                
                XYZ_i = cc.spd2XYZ(spec_i);
                XYZ_white = cc.spd2XYZ(spec_white);
                
                lab = CPR.ColorClass.XYZ2lab(XYZ_i,XYZ_white);
                
                obj.Lab(i,:) = lab;
            end
        end
        
        function spec_show (obj)
            %%SPEC_SHOW Show the 24 spectra in a figure

            % find max to size the plot
            ymax = max(max(obj.spec));
            
            for i = 1:24
                subplot(4,6,i)
                hold on
                plot(380:780,obj.spec(i,:))
                if 0
                axis([380 780 0 ymax])
                end
            end
        end

        function setXYZofTarget (obj,XYZ,target)
            % to replace setXYZ
            %assert(size(XYZ,1)==25)
            patch_n = target.patch_n;

            obj.XYZ = XYZ(1:patch_n,:);

            % decide reference white
            if size(XYZ,1) == target.patch_n + 1
                obj.XYZ_n = XYZ(patch_n+1,:);
            else 
                if size(XYZ,1) == target.patch_n
                    obj.XYZ_n = XYZ(target.white_patch,:);
                else
                    assert(1==0);
                end
            end

            for i = 1:patch_n
                obj.Lab(i,:) = CPR.ColorClass.XYZ2lab(obj.XYZ(i,:),obj.XYZ_n);
            end

            for i = 1:patch_n
                rgb = lab2rgb(obj.Lab(i,:));
                rgb = min(1,rgb);
                rgb = max(0,rgb);
                obj.rgb(i,:) = uint8(rgb*255);
            end

            obj.source = 'XYZ';
        end

        function setXYZ (obj,XYZ)
            % XYZ is 25x3
            %assert(size(XYZ,1)==25)

            obj.XYZ = XYZ(1:24,:);
            obj.XYZ_n = XYZ(25,:);

            for i = 1:24
                obj.Lab(i,:) = CPR.ColorClass.XYZ2lab(obj.XYZ(i,:),obj.XYZ_n);
            end

            for i = 1:24
                rgb = lab2rgb(obj.Lab(i,:));
                rgb = min(1,rgb);
                rgb = max(0,rgb);
                obj.rgb(i,:) = uint8(rgb*255);
            end

            obj.source = 'XYZ';
        end
 
       function setLAB (obj,LAB)
            % LAB is 24x3
            n = size(LAB,1);

            %assert(size(RGB,1)==24)

            obj.Lab = LAB(1:n,:);

            for i = 1:n
                obj.Lab(i,:) = LAB(i,:);
            end

            obj.source = 'LAB';
        end

         
        function setRGB (obj,RGB)
            % RGB is 24x3
            n = size(RGB,1);

            %assert(size(RGB,1)==24)

            obj.rgb = RGB(1:n,:);

            for i = 1:n
                obj.Lab(i,:) = rgb2lab(obj.rgb(i,:)/255,'colorspace','srgb','whitepoint','d65');
            end

            obj.source = 'RGB';
        end

        function calcRGB (obj)

            for i = 1:size(obj.Lab,1)
                rgb = lab2rgb(obj.Lab(i,:));
                rgb = min(1,rgb);
                rgb = max(0,rgb);
                obj.rgb(i,:) = uint8(rgb*255);
            end

        end

        function compare_with_ref (obj,ref,option)
            %%COMPARE_WITH_REF Compare this with the reference spectra
            % Use "option" to determine the scaling factor
            % none=1,0=auto, otherwise sc=option
            
            switch nargin
                case 2
                    r3 = 1;
                    
                case 3
                    if option == 0
                        % find ratio to align two spectra
                        r = obj.spec(:,50:end-50) ./ ref.spec(:,50:end-50);
                        r2 = mean(r,2);
                        r3 = mean(r2);
                    else
                        r3 = option;
                    end
            end
            
            sprintf('Adjusting reference by %fX',r3)
            
            % find max to size the plot
            ymax = max(max(obj.spec));
            
            for i = 1:24
                subplot(4,6,i)
                hold on
                plot(380:780,obj.spec(i,:),'r')
                plot(380:780,ref.spec(i,:)*r3,'b:')
                axis([380 780 0 ymax])
            end
        end
        
    end
end

