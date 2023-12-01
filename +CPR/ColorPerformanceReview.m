%% Color Performance Review (CPR) tool
% v2
% 12/1/2023

% called by cprapp with cpr:
%   one_pager()
%   K_number
%   MakeName
%   DeviceName
%   DUTName
%   input_colorspace
%   datapath

classdef ColorPerformanceReview < handle
    
    properties
        
        color_target          % test target id: 1, 2, 3
        color_target_name     % test target name: ColorChecker, ColorGauge, RezChecker
        all_patches           % list of all patches that are not empty
        grayscale_patches     % list of grayscale patches
        chromatic_patches     % list of chromatic patches
        all_patches_mask      % mask of all patches that are not empty
        grayscale_patches_mask % mask of grayscale patches
        chromatic_patches_mask % mask of chromatic patches
        
        devDataOriginal       % measurement data of the original (reference)
        devDataReproduced     % measurement data of the device (subject device)
        
        labOriginal       % LAB values of the original
        labReproduced     % LAB values of the device
        
        rgbOriginal       % RGB values of the original
        rgbReproduced     % RGB values of the device
        
        lchOriginal       % LCH values of the original
        lchReproduced     % LCH values of the device
        
        ilablchJoint      % joint table of id, L*, a*, b*, L*, c, h
        
        order_l           % order in lightness
        order_c           % order in chroma
        order_h           % order in hue
        order_a           % order in a*
        order_b           % order in b*
        
        im12_l6           % 1D order plot of the gray patches in lightness
        im12_l            % 1D order plot of the chromatic patches in lightness
        im12_c            % 1D order plot of the chromatic patches in chroma
        im12_h            % 1D order plot of the chromatic patches in hue
        
        kendall_l6        % Kendall tau-A of the gray patches in lightness
        kendall_l         % Kendall tau-A of the chromatic patches in lightness
        kendall_c         % Kendall tau-A of the chromatic patches in chroma
        kendall_h         % Kendall tau-A of the chromatic patches in hue
        
        target_type = 1;  % 24-patch, 30-patch
        patch_n           % number of patches, getting from CIELAB data
        
        CCE_name = [00 94 76]; % labels for different dE formulas
        
        dE_from_ref            % dE00, dE94, and dEab
        dEmean                 % stats for dE00
        dEstd                  % stats for dE00
        dEmax                  % stats for dE00
        dEmin                  % stats for dE00
        dEmedian               % stats for dE00
        
        CCEVectorOriginal      % all dE from truth as a vector
        CCEVectorReproduced    % all dE from device as a vector
        
        CCE00TableOriginal     % all dE00 from truth as a table for the original
        CCE94TableOriginal     % all dE94 from truth as a table for the original
        CCEabTableOriginal     % all dEab from truth as a table for the original
        
        CCE00TableReproduced   % all dE00 from truth as a table for the reproduced
        CCE94TableReproduced   % all dE94 from truth as a table for the reproduced
        CCEabTableReproduced   % all dEab from truth as a table for the reproduced
        
        CCE_ratio              % CCE ratios for dE00, dE94, and dEab
        CCE_mean               % mean CCE ratios for dE00, dE94, and dEab
        
        datapath               % path for retrieving sample endoscopic images
        
        K_number               % K number
        MakeName               % manufacturer name
        DeviceName             % device name
        DUTName                % device under test
        input_colorspace       % XYZ, LAB, RGB
        
    end
    
    %% Constructor
    methods
        
        function obj = ColorPerformanceReview (devDataRef, devDataSub)
            % Constructor of the ColorPerformanceReview class
            %   devDataRef: a DeviceData object for the reference dataset
            %   devDataSub: a DeviceData object for the device dataset
            
            [folder, name, ext] = fileparts(which('CPR.ColorPerformanceReview'));                 % where the data file is
            obj.datapath = [folder];                     % store the datapath
            
            switch size(devDataRef.Lab,1)
                case 24
                    obj.color_target = CPR.ColorCheckerClass;               % default ColorChecker
                case 30
                    obj.color_target = CPR.ColorGauge;               % default ColorChecker
                case 42
                    obj.color_target = CPR.RezChecker;               % default ColorChecker
                case 20
                    obj.color_target = CPR.SFRPlus;               % default ColorChecker
                otherwise
                    disp '!!!Color Target Unknown!!!'
                    return
            end
            
            patchid = [1:obj.color_target.patch_n]';
            
            obj.all_patches = obj.color_target.all_patches';
            obj.grayscale_patches = obj.color_target.grayscale_patches';
            obj.chromatic_patches = obj.color_target.chromatic_patches';
            
            obj.all_patches_mask = ismember(patchid,obj.all_patches);
            obj.grayscale_patches_mask = ismember(patchid,obj.grayscale_patches);
            obj.chromatic_patches_mask = ismember(patchid,obj.chromatic_patches);
            
            assert(size(devDataRef.Lab,1)==size(devDataSub.Lab,1)); % check whether Ref and Sub have the same size
            
            obj.patch_n = size(devDataRef.Lab,1);                   % save the patch count
            obj.input_colorspace = devDataRef.source;
            obj.color_target_name = obj.color_target.name;
            
            assert(obj.patch_n == obj.color_target.patch_n);        % check whether Ref and ColorChecker have the same patch count
            
            % save the original DeviceData
            obj.devDataOriginal = devDataRef;
            obj.devDataReproduced = devDataSub;
            
            % save the LAB
            obj.labOriginal = devDataRef.Lab;
            obj.labReproduced = devDataSub.Lab;
            
            % calculate the RGB
            obj.rgbOriginal = lab2rgb(obj.labOriginal)*255;
            obj.rgbReproduced = lab2rgb(obj.labReproduced)*255;
            
            obj.rgbOriginal = min(obj.rgbOriginal,255);
            obj.rgbOriginal = max(obj.rgbOriginal,0);
            
            obj.rgbReproduced = min(obj.rgbReproduced,255);
            obj.rgbReproduced = max(obj.rgbReproduced,0);
            
            % calculate the LCH
            obj.lchOriginal = obj.lab2lch(obj.labOriginal);
            obj.lchReproduced = obj.lab2lch(obj.labReproduced);
            
            obj.sort_labch;
            
            obj.calcualte_concordance;
            
            obj.calculate_dE_from_truth;
            
            obj.calculate_CCE;
            
            obj.calculate_visual_1x24;
        end
    end

    %% Helper functions start with "calculate_"
    methods (Access = private)

        function sort_labch (obj)
            % Sort patches by L*, a*, b*, c, and h
            
            patch_n = obj.patch_n;
            
            % construct ilabch
            ilablchOriginal = [[1:patch_n]' obj.labOriginal obj.lchOriginal];
            ilablchReproduced = [[1:patch_n]' obj.labReproduced obj.lchReproduced];
            
            ilablchJoint = [ilablchOriginal ilablchReproduced];
            
            obj.ilablchJoint = ilablchJoint;
            
            obj.order_l = helper(ilablchJoint,2,9);
            obj.order_a = helper(ilablchJoint,3,10);
            obj.order_b = helper(ilablchJoint,4,11);
            obj.order_c = helper(ilablchJoint,6,13);
            obj.order_h = helper(ilablchJoint,7,14);
            
            return
            
            function order = helper (tab,col_ref,col_subject)
                % Sort "tab" by col_ref and by col_subject
                
                sorted_for_ref = sortrows(tab,col_ref);
                order(:,1) = sorted_for_ref(:,1);
                
                sorted_for_subject = sortrows(tab,col_subject);
                order(:,2) = sorted_for_subject(:,8);
            end
        end
        
        function calcualte_concordance (obj)
            % Calculate kendall_l6, kendall_l, kendall_c, and kendall_h
            
            x = obj.labOriginal(obj.grayscale_patches_mask,1);
            y = obj.labReproduced(obj.grayscale_patches_mask,1);
            obj.kendall_l6 = obj.kendall_tau_a(x,y);
            
            x = obj.labOriginal(obj.chromatic_patches_mask,1);
            y = obj.labReproduced(obj.chromatic_patches_mask,1);
            obj.kendall_l = obj.kendall_tau_a(x,y);
            
            x = obj.lchOriginal(obj.chromatic_patches_mask,2);
            y = obj.lchReproduced(obj.chromatic_patches_mask,2);
            obj.kendall_c = obj.kendall_tau_a(x,y);
            
            x = obj.lchOriginal(obj.chromatic_patches_mask,3);
            y = obj.lchReproduced(obj.chromatic_patches_mask,3);
            obj.kendall_h = obj.kendall_tau_a(x,y);
            
        end
        
        function calculate_dE_from_truth (obj)
            % Calculate dE from truth
            
            dE = zeros(obj.patch_n,3);
            
            % mask = [1:obj.white_patch-1 obj.white_patch+1:obj.patch_n];  % remove #19 white??
            mask = obj.all_patches;
            
            for i = 1:obj.patch_n
                
                [dE00 dE94 dEab] = CPR.ColorClass.LAB2dE(obj.labOriginal(i,:)',obj.labReproduced(i,:)');
                
                dE(i,1:3) = [dE00 dE94 dEab];
                
            end
            
            obj.dE_from_ref = dE;
            
            obj.dEmean = mean(dE(mask));
            obj.dEstd = std(dE(mask));
            obj.dEmax = max(dE(mask));
            obj.dEmin = min(dE(mask));
            obj.dEmedian = median(dE(mask));
        end
        
        function calculate_CCE (obj)
            % Analyze CCE
            
            devOriginal = obj.devDataOriginal;
            devReproduced = obj.devDataReproduced;
            
            n = numel(obj.all_patches);
            
            truth(1:n) = CPR.ColorClass([0 0 0]);
            measured(1:n) = CPR.ColorClass([0 0 0]);
            
            for i = 1:n
                ii = obj.all_patches(i);
                truth(i) = devOriginal.Lab(ii,:);
                measured(i) = devReproduced.Lab(ii,:);
            end
            
            % dE pairwise
            obj.CCEVectorOriginal = zeros(n*(n-1)/2,3);
            obj.CCEVectorReproduced = zeros(n*(n-1)/2,3);
            
            obj.CCE00TableOriginal = zeros(n,n);
            obj.CCE94TableOriginal = zeros(n,n);
            obj.CCEabTableOriginal = zeros(n,n);
            
            obj.CCE00TableReproduced = zeros(n,n);
            obj.CCE94TableReproduced = zeros(n,n);
            obj.CCEabTableReproduced = zeros(n,n);
            
            % the color patches considered for evaluation
            k = 0;
            for ii = 1:n
                for jj = ii+1:n
                    k = k + 1;
                    i = ii;
                    j = jj;
                    
                    % do reference
                    [dE00 dE94 dEab] = CPR.ColorClass.LAB2dE (truth(i).lab', truth(j).lab');
                    obj.CCEVectorOriginal(k,:) = [dE00 dE94 dEab]';
                    
                    obj.CCE00TableOriginal(i,j) = dE00;
                    obj.CCE94TableOriginal(i,j) = dE94;
                    obj.CCEabTableOriginal(i,j) = dEab;
                    
                    % do device
                    [dE00 dE94 dEab] = CPR.ColorClass.LAB2dE (measured(i).lab', measured(j).lab');
                    obj.CCEVectorReproduced(k,:) = [dE00 dE94 dEab]';
                    
                    obj.CCE00TableReproduced(i,j) = dE00;
                    obj.CCE94TableReproduced(i,j) = dE94;
                    obj.CCEabTableReproduced(i,j) = dEab;
                end
            end
            
            obj.CCE_ratio = obj.CCEVectorReproduced ./ obj.CCEVectorOriginal;
            obj.CCE_mean = mean(obj.CCE_ratio);
        end
        
        function calculate_visual_1x24 (obj)
            % Order in lightness, chroma, and hue - 1D view
            
            % 'Lightness6'
            obj.im12_l6 = helper(obj.order_l,obj.grayscale_patches,obj.kendall_l6);
            
            % 'Lightness'
            obj.im12_l = helper(obj.order_l,obj.chromatic_patches,obj.kendall_l);
            
            % 'Chroma'
            obj.im12_c = helper(obj.order_c,obj.chromatic_patches,obj.kendall_c);
            
            % 'Hue'
            obj.im12_h = helper(obj.order_h,obj.chromatic_patches,obj.kendall_h);
            
            return
            
            function im12 = helper (order,range,tau)
                cc1 = CPR.ColorCheckerClass_1x24;
                
                mask = ismember(order(:,1),range);
                idx1 = order(mask,1);
                rgb_l24 = obj.rgbOriginal(idx1,:);
                im1 = cc1.getImage_with_index(rgb_l24,idx1);
                
                mask = ismember(order(:,2),range);
                idx2 = order(mask,2);
                rgb_l24 = obj.rgbReproduced(idx2,:);
                im2 = cc1.getImage_with_index(rgb_l24,idx2);
                
                im12 = stitching_two_colorbars(im1,im2,idx1,idx2);
                
                %
                % add tau
                %
                tau_str = sprintf('Tau %.2f',tau);
                im12 = insertText(im12,[10 size(im12,1)/2],tau_str,TextColor='black',FontSize=24,BoxColor='green');
                
                return
             
            end
            
            function im12 = stitching_two_colorbars (im1,im2,idx1,idx2)
                % Stitch two colorbars together with lines

                im3 = im1;

                im3 = im3*0.01;

                im12 = [im1 ; im3; im2];

                n = length(idx1);

                idx2two = [idx2 [1:n]'];
                idx2twosorted = sortrows(idx2two,1);

                imw = size(im12,2);
                imh = size(im12,1);

                boxw = imw/n;

                posx = [1:length(idx1)]*boxw - boxw/2 ;
                posy = [boxw*1 boxw*2];

                for location1 = 1:n
                    patch_id = idx1(location1);
                    location2 = find(idx2 == patch_id);

                    im12 = insertShape(im12,'Line',[posx(location1) posy(1) posx(location2) posy(2)],'Color','Yellow','LineWidth',5);
                end

                % for debugging
                %             clf
                %             imshow(im12)

                return
            end

        end
       
    end
    
    %% Functions used by this class
    %    keep them public so that others can use
    methods
        
        function tau_a = kendall_tau_a (obj,x,y)
            % Kendall rank correlation coefficient
            % https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient
            
            n = length(x);
            pairs = 0;
            concordant = 0;
            discordant = 0;
            for i = 1:n
                for j = i+1:n
                    pairs = pairs + 1;
                    if sign(x(i)-x(j)) * sign(y(i)-y(j)) >= 0
                        concordant = concordant+1;
                    else
                        discordant = discordant+1;
                    end
                end
            end
            %[concordant discordant pairs]
            tau_a = (concordant-discordant)/pairs;
        end
        
        function [p, yfit, rsq] = linear_regression (obj,x,y)
            % linear regression of a scatter plot
            
            p = polyfit(x,y,1);
            yfit = polyval(p,x);
            yresid = y - yfit;
            SSresid = sum(yresid.^2);
            SStotal = (length(y)-1) * var(y);
            rsq = 1 - SSresid/SStotal;
        end
        
        function lch = lab2lch (obj,lab)
            % Convert CIELAB to CIELCH
            
            l = lab(:,1);
            a = lab(:,2);
            b = lab(:,3);
            lch(:,1) = l;
            lch(:,2) = (a.^2 + b.^2).^0.5;
            lch(:,3) = atan2d(b,a);
        end
        
        function save_figure (obj,resizeto,filename)
            % Save the plot with the given size
            
            my_gcf = gcf;
            my_gcf.Units = 'pixels';
            my_gcf.Position = resizeto;
            saveas(gcf,['fig/' filename])
        end
        
    end

    %% one-pager functions
    %
    methods (Access = private)

        function [im0, im1, im2] = evaluate_visual_XYZ_im (obj)
            % Visualize CIEXYZ data
            % Figure (a1)
            %   If the reproduced data were provided in the CIEXYZ color space (i.e., tristimulus),
            %   an optional Figure 8(a1) shows the color target by using CIE D65 as the reference white.
            %   The purpose is to check any excessive color shift caused by the endoscope light source and/or the device.            

            XYZ_d65 = rgb2xyz([1 1 1]);                                    % take the white point (D65)
            Y_d65 = XYZ_d65(2);                                            % retrieve the luminance Y

            XYZ1_orig = obj.devDataOriginal.XYZ;                           % 1: original
            Y1_orig = XYZ1_orig(19,2);                                     %   take the white patch #19 from ColorChecker
            XYZ1 = XYZ1_orig / Y1_orig * Y_d65;                            %   normalize

            XYZ2_orig = obj.devDataReproduced.XYZ;                         % 2: reproduced
            Y2_orig = XYZ2_orig(19,2);                                     %   take the white patch #19 from ColorChecker
            XYZ2 = XYZ2_orig / Y2_orig * Y_d65;                            %   normalize

            %
            %
            %
            ccc = CPR.ColorCheckerClass;

            im0 = ccc.getImage_with_index;                                       % create a colorchecker with the default RGB values

            rgb1 = xyz2rgb(XYZ1,'ColorSpace','srgb','WhitePoint','d65')*255;     % create a colorchecker with XYZ1
            im1 = ccc.getImage(rgb1);

            rgb2 = xyz2rgb(XYZ2,'ColorSpace','srgb','WhitePoint','d65')*255;     % create a colorchecker with XYZ2
            im2 = ccc.getImage(rgb2);
        end

        function evaluate_visual (obj, ha0, ha1, ha2)
            % Visualize CIEXYZ data
            % Figure (a2)(a3)(a4)
            %   ha0: original
            %   ha1: reference
            %   ha2: reproduced
            
            [im0, im1, im2] = evaluate_visual_im;
            
            axes(ha0)
            imshow(im0)
            title('Original ColorChecker')
            
            axes(ha1)
            imshow(im1)
            title('Reference Data, LAB')
            
            axes(ha2)
            imshow(im2)
            title('Reproduced Data, LAB')

            return 

            function [im0, im1, im2] = evaluate_visual_im
                % Evaluate visual im

                ccc = obj.color_target;

                im0 = ccc.getImage_with_index;
                im1 = ccc.getImage(obj.rgbOriginal);
                im2 = ccc.getImage(obj.rgbReproduced);
            end

        end
        
        function im2 = simulate_im (obj,im1)
            % Interpolate the image to be simulated
            % Figure (a5)(a6)
            %   In addition to the color target, the user can provide an endoscopic image (Figure 8(a5)) for the tool to predict the device's output (Figure 8(a6)). 
            %   The prediction is based on linear interpolation of the patches. The purpose is to provide a quick visual assessment of the color reproduction of a real scene.
            
            rgb1 = reshape(im1,size(im1,1)*size(im1,2),3);
            lab1 = rgb2lab(rgb1);
            
            labref = obj.labOriginal;
            labtarget = obj.labReproduced;
            
            Fl = scatteredInterpolant(labref,labtarget(:,1));
            Fa = scatteredInterpolant(labref,labtarget(:,2));
            Fb = scatteredInterpolant(labref,labtarget(:,3));
            
            mymethod = 'linear';
            
            Fl.Method = mymethod;
            Fa.Method = mymethod;
            Fb.Method = mymethod;
            
            lab2 = [Fl(lab1(:,1),lab1(:,2),lab1(:,3)) Fa(lab1(:,1),lab1(:,2),lab1(:,3)) Fb(lab1(:,1),lab1(:,2),lab1(:,3))];
            rgb2 = lab2rgb(lab2);

            im2 = reshape(rgb2,size(im1,1),size(im1,2),3);
        end
        
        function evaluate_dE_from_truth (obj,ha1,ha2)
            % Evaluate dE from the ground truth
            % Figure (a7)(a8)
            
            metric = 1;
            
            axes(ha1)                                                      % Figure (a7)    
            h1 = bar(obj.all_patches,obj.dE_from_ref(obj.all_patches_mask,metric),'FaceColor','flat');
            h1.CData(:,:) = obj.rgbOriginal(obj.all_patches_mask,:)/255;
            xlabel('Patch #')
            ylabel('Color Difference \DeltaE_{00}')
            title(sprintf('\\mu=%0.2f, \\delta=%0.2f', ...
                obj.dEmean(metric), ...
                obj.dEstd(metric)))
            
            %
            % problem: xlabel in bar chart was clipped
            % solution: https://www.mathworks.com/matlabcentral/answers/37422-ensure-x-label-never-cut-off
            %
            
            op = get(gca,'OuterPosition');
            set(gca,'OuterPosition',op+[0 0.1 0 0])
            
            axes(ha2)                                                      % Figure (a8)
            set(gca,'TickLabelInterpreter','latex');
            
            %            boxplot(obj.dE_from_ref(obj.all_patches_mask,metric))
            boxplot(obj.dE_from_ref(obj.all_patches_mask,1:3),{'dE00','dE94','dE76'});
            
            ylabel('Color Difference \DeltaE_{00}')
            
            % to align the height between two subplots
            %            axis([0 4 h1.Parent.YLim(1) h1.Parent.YLim(2)])
            
            op = get(gca,'OuterPosition');
            set(gca,'OuterPosition',op+[0 0.1 0 0])
            
            title(sprintf('min=%0.2f, median=%0.2f, max=%0.2f', ...
                obj.dEmin(metric), ...
                obj.dEmedian(metric), ...
                obj.dEmax(metric)))
        end        

        function evaluate_linearity_lightness_6 (obj, ha)
            % show linearity plot for lightness 6 patches
            % Figure (b1)
            
            axes(ha);
            
            l_original = obj.labOriginal(obj.grayscale_patches_mask,1);
            l_reproduced = obj.labReproduced(obj.grayscale_patches_mask,1);
            id = obj.grayscale_patches;
            
            [p, yfit, rsq] = obj.linear_regression(l_original,l_reproduced);
            
            obj.linearity_plot(l_original,l_reproduced,id,[0 100 0 100],p);
            
            xlabel('L* reference')
            ylabel('L* measured')
            title('Lightness Order')            
        end
        
        function evaluate_linearity_lightness (obj, ha)
            % show linearity plot for lightness
            % Figure (b2)
            
            axes(ha);
            
            l_original = obj.labOriginal(obj.chromatic_patches_mask,1);
            l_reproduced = obj.labReproduced(obj.chromatic_patches_mask,1);
            id = obj.chromatic_patches;
            
            [p, yfit, rsq] = obj.linear_regression(l_original,l_reproduced);
            
            obj.linearity_plot(l_original,l_reproduced,id,[0 100 0 100],p);
            
            xlabel('L* reference')
            ylabel('L* measured')
            title('Lightness Order')            
        end
        
        function linearity_plot (obj, x, y, id, box, p)
            % show one linearity plot
            %   a helper function used by Figure (b3)(b4)
            
            hold on
            for i=1:length(id)
                rgb256 = obj.rgbOriginal(id(i),:);
                rgb = rgb256 / 255.0;
                rgb = min(1,rgb);
                rgb = max(0,rgb);
                plot(x(i),y(i),'o','MarkerFaceColor',rgb,'MarkerEdgeColor',rgb)
            end
            text(x,y+1,num2str(id))
            
            axis(box)
            plot(box(1:2),box(3:4),'-k')   % add the identity line
            yfit = polyval(p,box(1):box(2));
            plot(box(1):box(2),yfit,':r')   % add the fitting line
            
            fitting_eq = sprintf('y=%.2fx+%.2f',p(1),p(2));
            text(box(1)+(box(2)-box(1))/10,box(3)+(box(4)-box(3))/10*9,fitting_eq)
            
            axis square
            grid on
        end

        function evaluate_linearity_hue (obj,ha)
            % show linearity plot for hue
            % Figure (b3)
            
            axes(ha);
            
            h_original = obj.lchOriginal(obj.chromatic_patches_mask,3);
            h_reproduced = obj.lchReproduced(obj.chromatic_patches_mask,3);
            id = obj.chromatic_patches;
            
            [p, yfit, rsq] = obj.linear_regression(h_original,h_reproduced);
            
            obj.kendall_h = obj.kendall_tau_a(h_original,h_reproduced);
            
            obj.linearity_plot(h_original,h_reproduced,id,[-180 180 -180 180],p);
            
            xlabel('Hue reference')
            ylabel('Hue measured')
            title('Hue Order')           
        end

        function evaluate_linearity_chroma (obj,ha)
            % show linearity plot for chroma
            % Figure (b4)
            
            axes(ha);
            
            c_original = obj.lchOriginal(obj.chromatic_patches_mask,2);
            c_reproduced = obj.lchReproduced(obj.chromatic_patches_mask,2);
            id = obj.chromatic_patches;
            
            [p, yfit, rsq] = obj.linear_regression(c_original,c_reproduced);
            
            obj.kendall_c = obj.kendall_tau_a(c_original,c_reproduced);
            
            obj.linearity_plot(c_original,c_reproduced,id,[0 100 0 100],p);
            
            xlabel('Chroma reference')
            ylabel('Chroma measured')
            title('Chroma Order')            
        end
        
        function order_1D_l6 (obj,ha)
            % Order 1D lightness 6
            % Figure (b5)
            
            axes(ha)
            obj.layout_1D(obj.order_l, obj.grayscale_patches, obj.kendall_l6);
        end
        
        function order_1D_l (obj,ha)
            % Order 1D lightness
            % Figure (b6)
            
            axes(ha)
            obj.layout_1D(obj.order_l, obj.chromatic_patches, obj.kendall_l);
        end
        
        function order_1D_c (obj,ha)
            % Order 1D 'Chroma'
            % Figure (b7)
            
            axes(ha)
            obj.layout_1D(obj.order_c, obj.chromatic_patches, obj.kendall_c);
        end
        
        function order_1D_h (obj,ha)
            % Order 1D 'Hue'
            % Figure (b8)
            
            axes(ha)
            obj.layout_1D(obj.order_h, obj.chromatic_patches, obj.kendall_h);
        end
        
        function layout_1D (obj,order,range,tau)
            % Order 1D
            % Layout function for Figures (b5)(b6)(b7)(b8)
            
            cc1 = CPR.ColorCheckerClass_1x24;
            
            mask = ismember(order(:,1),range);
            idx1 = order(mask,1);
            rgb1 = obj.rgbOriginal(idx1,:);
            
            mask = ismember(order(:,2),range);
            idx2 = order(mask,2);
            rgb2 = obj.rgbReproduced(idx2,:);
            
            offset = 0.4;
            row_from = 0.5+offset;
            row_to = 0.5-offset;
            markersize = 15;
            col_gap = 1/20;
            
            %hf = gca;
            %hf.Units = 'Inches';
            
            hold on
            
            for i = 1:numel(idx1)
                mycolor = rgb1(i,:)/255;
                mystr = sprintf('%02d',idx1(i));
                
                xsize = 1/21;
                ysize = 1/5;
                pos = [i*col_gap row_from xsize ysize];
                
                %text(pos(1),pos(2),mystr);
                %upper{i} = annotation('rectangle',pos,'FaceColor',mycolor);
                plot(pos(1),pos(2),'s','MarkerSize',markersize,'MarkerFaceColor',mycolor,'MarkerEdgeColor',mycolor);
                upper{i} = pos;
            end
            
            for i = 1:numel(idx2)
                mycolor = rgb2(i,:)/255;
                mystr = sprintf('%02d',idx2(i));
                pos = [i*col_gap row_to xsize ysize];
                %                 ha = annotation('textbox',pos, ...
                %                     'String',mystr,...
                %                     'BackgroundColor',mycolor);
                % lower{i} = annotation('rectangle',pos,'FaceColor',mycolor);
                plot(pos(1),pos(2),'s','MarkerSize',markersize,'MarkerFaceColor',mycolor,'MarkerEdgeColor',mycolor);
                lower{i} = pos;
            end
            
            for location1 = 1:numel(idx1)
                i_from = idx1(location1);
                i_to = idx2(location1);
                
                patch_id = idx1(location1);
                location2 = find(idx2 == patch_id);
                mycolor = rgb1(location1,:)/255;
                
                %                 x1 = upper{location1}.Position(1) + upper{location1}.Position(3)/2 ;
                %                 y1 = upper{location1}.Position(2);
                %                 x2 = lower{location2}.Position(1) + upper{location1}.Position(3)/2;
                %                 y2 = lower{location2}.Position(2) + ysize;
                
                x1 = upper{location1}(1);
                y1 = upper{location1}(2);
                x2 = lower{location2}(1);
                y2 = lower{location2}(2);
                plot([x1 x2],[y1 y2],'-','Color',mycolor,'LineWidth',2);
                
                % ha = annotation('line',[x1 x2],[y1 y2]);
                %ha.Parent = hf.CurrentAxes;
                
            end
            
            axis off
            axis([0 1 0 1.2])
            %im12 = obj.stitching_two_colorbars(im1,im2,idx1,idx2);
            
            %
            % add tau
            %
            %tau_str = sprintf('Tau %.2f',tau);
            %im12 = insertText(im12,[10 size(im12,1)/2],tau_str,TextColor='black',FontSize=24,BoxColor='green');
            
        end
         
        function show3dquiver_select (obj,mask,ha)
            % Show color transfer of a subset in CIELAB
            % Figure (c1)(c3)(c4)(c5)
            
            axes(ha)                                                       % select the axes
            
            lab1 = obj.labOriginal(mask,:);                                % convert LAB to RGB for painting
            lab2 = obj.labReproduced(mask,:);
            
            rgb1 = lab2rgb(lab1);
            rgb1 = min(1,rgb1);
            rgb1 = max(0,rgb1);
            
            rgb2 = lab2rgb(lab2);
            rgb2 = min(1,rgb2);
            rgb2 = max(0,rgb2);
            
            dif = lab2-lab1;                                               % convert LAB to RGB for painting
            
            hold on
            
            quiver3(lab1(:,2),lab1(:,3),lab1(:,1),dif(:,2),dif(:,3),dif(:,1),0,'k.')            % draw quivers
            
            for i = 1:length(mask)                                                              % add end points  
                plot3(lab1(i,2),lab1(i,3),lab1(i,1),'o','MarkerFaceColor',rgb1(i,:),'MarkerEdgeColor',rgb1(i,:))
                plot3(lab2(i,2),lab2(i,3),lab2(i,1),'x','MarkerFaceColor',rgb2(i,:),'MarkerEdgeColor',rgb2(i,:))
            end
            
            zlabel('L*')
            xlabel('a*')
            ylabel('b*')
            grid on
            axis equal
        end

        function show2dquiver_polar_select (obj,mask,ha)
            % Show hue/chroma transfer of a subset in CIELAB on a polar chart
            % Figure (c2)(c6)
            
            axes(ha)                                                       % select the axes
            
            show3dquiver_select (obj,mask,ha)                              % show the quivers
            
            line_color = 0.8;
            
            plot_circle(25)                                                % draw 4 concentric rings
            plot_circle(50)
            plot_circle(75)
            plot_circle(100)
            
            for i = 0:11                                                   % draw 12 radii
                th = 2*pi/12 * i;
                plot_radius(th,100)
            end
            
            grid off                                                       % turn off grid 
                        
            return
            
            function plot_circle (r)
                % plot a ring
                
                th = 0:0.01:2*pi;
                plot(r*cos(th),r*sin(th),'-','Color',[1 1 1]*line_color)
            end
            
            function plot_radius (th,r)
                % plot a radius of length r at angle th
                
                plot([0 r*cos(th)],[0 r*sin(th)],'-','Color',[1 1 1]*line_color)
            end
        end
        
        function evaluate_CCE (obj,metric,ha)
            % Show dE_in vs. dE_out
            %   metric: 1=dE00, 2=dE94, 3=dEab
            % Figure (c7)
            
            axes(ha)
            
            x = obj.CCEVectorOriginal(:,metric);
            y = obj.CCEVectorReproduced(:,metric);
            [p,yfit,rsq] = obj.linear_regression(x,y);
            
            yoverx = y./x;
            cce_good = nnz(yoverx >= 1);
            cce_all = numel(yoverx);
            cce_good_percent = cce_good / cce_all * 100;
            
            box = [0 100 0 100];
            
            hold on
            
            switch metric
                case 1
                    CCETableIn = obj.CCE00TableOriginal;
                    CCETableOut = obj.CCE00TableReproduced;
                case 2
                    CCETableIn = obj.CCE94TableOriginal;
                    CCETableOut = obj.CCE94TableReproduced;
                case 3
                    CCETableIn = obj.CCEabTableOriginal;
                    CCETableOut = obj.CCEabTableReproduced;
            end
            
            n = numel(obj.all_patches);
            for i = 1:n
                for j = i+1:n
                    dein = CCETableIn(i,j);
                    deout = CCETableOut(i,j);
                    colorin = obj.rgbOriginal(i,:)/255;
                    colorout = obj.rgbReproduced(i,:)/255;
                    plot(dein,deout,'_','MarkerEdgeColor',colorin);
                    plot(dein,deout,'|','MarkerEdgeColor',colorout);
                end
            end
            
            % plot(x,y,'o')
            
            plot(box(1:2),box(3:4),'--k')   % add the identity line
            
            %
            % add fitting data
            %
            if 0
                yfit = polyval(p,box(1):box(2));
                plot(box(1):box(2),yfit,':r')   % add the fitting line
                
                fitting_eq = sprintf('y=%.2fx+%.2f',p(1),p(2));
                text(box(1)+(box(2)-box(1))/10,box(3)+(box(4)-box(3))/10*9,fitting_eq);
            end
            
            %
            % add good CCE percent
            %
            if 1
                cce_eq = sprintf('%.2f%%',cce_good_percent);
                text(box(1)+(box(2)-box(1))/10,box(3)+(box(4)-box(3))/10*9,cce_eq);
            end
            
            grid on
            axis square
            axis(box)
            xlabel('Input \Delta{E}')
            ylabel('Output \Delta{E}')
            title(sprintf('CCE_{%02d}=%.2f%c%0.2f', ...
                obj.CCE_name(metric), ...
                obj.CCE_mean(metric), ...
                char(177),...
                std(obj.CCE_ratio(:,metric))))
            %legend('Measured','Identity')
            %legend('Location','Southeast')

        end

    end
    
    %% Funcitons for CPR GUI
    methods
        
        %% one_pager
        % <<fig 8 DeviceA_one_pager_annotated.png>>

        function one_pager (cpr, fnout)
            % Generate the one-page report
            %   fnout: output filename
            
            if nargin < 2                                                  % if the output filename is not supplied
                fnout = 'one_pager.png';                                   %   use the default output filename
            end
            
            clf                                                            % start the canvas
            hg = gcf;
            set(hg,'Visible','on')
            
            td = tiledlayout(8,6,'TileSpacing','tight');                   % divide the canvas into 8x6 cells
            
            title(td,sprintf('%s, %s, %s, %s, %d, %s, %s', ...             % set the title
                cpr.K_number, ...
                cpr.MakeName, ...
                cpr.DeviceName, ...
                cpr.color_target_name, ...
                cpr.patch_n, ...
                cpr.input_colorspace, ...
                cpr.DUTName ...
                ))
            
            
            %% Visual verification and absolute color error
            %
            %

            % If the reproduced data were provided in the CIEXYZ color space (i.e., tristimulus), an optional Figure (a1) shows the color target by using CIE D65 as the reference white. 
            % The purpose is to check any excessive color shift caused by the endoscope light source and/or the device.
            nexttile(1,[2 1]);
            if strcmp(cpr.devDataOriginal.source,'XYZ')                    % show it only when the input data are XYZ
                                                                           % ! This part will never be reached from CPR GUI because CPR GUI passes CIELAB data. ! 
                [im0, im1, im2] = cpr.evaluate_visual_XYZ_im;              
                imshow(im2)
                title('Reproduced Data, XYZ')
            else
                axis off
            end
            
            % Figure (a2) shows the original color target with each patch numbered. 
            ha0 = nexttile(2,[2 1]);                                       % 2x1 cell at [2 8]

            % Figure (a3) and (a4) show the color target based on the converted CIELAB data for the reference and reproduced data, respectively. 
            ha1 = nexttile(13,[2 1]);                                      % 2x1 cell at [13 19]
            ha2 = nexttile(14,[2 1]);                                      % 2x1 cell at [14 20]
            
            cpr.evaluate_visual(ha0,ha1,ha2);
            
            %% Endoscopic image simulation
            % the user can provide an endoscopic image (Figure (a5)) for the tool to predict the device's output (Figure (a6))
            
            %            filename = [obj.datapath '/sample_polyp.png'];
            filename = ['sample_polyp.png'];
            
            im1 = imread(filename);
            im2 = cpr.simulate_im (im1);
            
            % Fig a5
            ha1 = nexttile(25,[2 1]);
            axes(ha1);
            imshow(im1);
            title('Original Scene')

            % Fig a6
            ha2 = nexttile(26,[2 1]);
            axes(ha2);
            imshow(im2);
            title('Predicted Device Output')
            
            %% Absolute dE
            %

            % Figure (a7) is a bar chart showing the absolute color errors between the device output and the ground truth for each patc
            ha1 = nexttile(37,[2 1]);

            % Figure (a8) is boxplot showing the color differences calculated by using the ∆E_00, ∆E_94, and ∆E_76 formulas.
            ha2 = nexttile(38,[2 1]);

            cpr.evaluate_dE_from_truth(ha1,ha2);
            
            %% Preservation of the patch order in lightness, hue, and chroma
            %
            %

            %% Order 2D (linearity)
            %    evaluating the linearity 

            % Figure (b1) shows lightness order of the gray patches
            ha1 = nexttile(3,[2 1]);
            cpr.evaluate_linearity_lightness_6(ha1); title(sprintf('Lightness order of gray patches'))
            
            % Figure (b2) shows lightness order of the chromatic patches
            ha2 = nexttile(4,[2 1]);
            cpr.evaluate_linearity_lightness(ha2); title(sprintf('Lightness order of chromatic patches'))
            
            % Figure (b3) shows hue order of the chromatic patches
            ha3 = nexttile(15,[2 1]);
            cpr.evaluate_linearity_hue(ha3); title(sprintf('Hue order of chromatic patches'))
            
            % Figure (b4) shows chroma order of the chromatic patches
            ha4 = nexttile(16,[2 1]);
            cpr.evaluate_linearity_chroma(ha4); title(sprintf('Chroma order of chromatic patches'))
            
            %% Order 1D
            %    evaluating the order 
            
            % Figure (b5) shows lightness order of the gray patches
            ha1 = nexttile(27,[1 2]); cpr.order_1D_l6(ha1); title(sprintf('Lightness order of gray patches,{\\tau}=%.2f',cpr.kendall_l6))

            % Figure (b6) shows lightness order of the chromatic patches
            ha2 = nexttile(33,[1 2]); cpr.order_1D_l(ha2); title(sprintf('Lightness order of chromatic patches,{\\tau}=%.2f',cpr.kendall_l))

            % Figure (b7) shows hue order of the chromatic patches
            ha3 = nexttile(39,[1 2]); cpr.order_1D_h(ha3); title(sprintf('Hue order of chromatic patches,{\\tau}=%.2f',cpr.kendall_h))

            % Figure (b8) shows chroma order of the chromatic patches
            ha4 = nexttile(45,[1 2]); cpr.order_1D_c(ha4); title(sprintf('Chroma order of chromatic patches,{\\tau}=%.2f',cpr.kendall_c))
            
            %% 3D quiver
            % Figure (c1)-(c6) allow the user to visually examine the three-dimensional color transfer as vectors in the CIELAB color space. 

            % Figure (c1)
            ha1 = nexttile(5,[2 1]); cpr.show3dquiver_select(cpr.chromatic_patches,ha1); view([12 32]); axis equal; title('Color transfer')

            % Figure (c2)
            ha2 = nexttile(6,[2 1]); cpr.show2dquiver_polar_select(cpr.chromatic_patches,ha2); view([0 90]); axis([-100 100 -100 100 0 100]); title('Hue Chroma')

            % Figure (c3)
            ha3 = nexttile(17,[2 1]); cpr.show3dquiver_select(cpr.chromatic_patches,ha3); view([0 0]); axis square; axis([-100 100 -100 100 0 100]); title('L* a* plane')

            % Figure (c4)
            ha4 = nexttile(18,[2 1]); cpr.show3dquiver_select(cpr.chromatic_patches,ha4); view([90 0]); axis square; axis([-100 100 -100 100 0 100]); title('L* b* plane')

            % Figure (c5)
            ha5 = nexttile(29,[2 1]); cpr.show3dquiver_select(cpr.grayscale_patches,ha5); view([0 0]); axis square; axis([-100 100 -100 100 0 100]); title('Lightness gray')

            % Figure (c6)
            ha6 = nexttile(30,[2 1]); cpr.show2dquiver_polar_select(cpr.grayscale_patches,ha6); view([0 90]); axis square; axis([-30 30 -30 30 0 100]); title('Chromaticity gray')
            
            %% CCE scatter
            % Figure (c7) is a scatter plot showing the relationship between the input ∆E, output ∆E, and CCE
            ha1 = nexttile(41,[2 1]); cpr.evaluate_CCE(1,ha1);
            % ha2 = nexttile(24); cpr.evaluate_CCE(2,ha2);
            % ha3 = nexttile(27); cpr.evaluate_CCE(3,ha3);
            
            %% CCE boxplot
            % Figure (c8) is a box plot showing the CCE distributions using the ∆E_00, ∆E_94, and ∆E_76 formulas. 
            ha1 = nexttile(42,[2 1]);
            boxplot(cpr.CCE_ratio,{'CCE00','CCE94','CCE76'})
            ylabel('CCE (%)')
            title('CCE by different \DeltaE')
            
            %% Save the image
            %
            set(hg,'Position',[137 162 1714 939])
            saveas(hg,fnout)
            
        end
        
        function cpr = one_pager_from_LAB (LAB_r, LAB_s, K_number, MakeName, DeviceName, DUTName, fnout)
            % Generate one-pager from LAB
            
            LAB_Reference = CPR.DeviceData;
            LAB_Reference.setLAB(LAB_r);
            
            LAB_Subject = CPR.DeviceData;
            LAB_Subject.setLAB(LAB_s);
            
            cpr = CPR.ColorPerformanceReview(LAB_Reference,LAB_Subject);
            cpr.K_number = K_number;
            cpr.MakeName = MakeName;
            cpr.DeviceName = DeviceName;
            cpr.DUTName = DUTName;
            
            cpr.one_pager
            
            movefile('one_pager.png',fnout)
            
        end
                
        function cpr = one_pager_from_XYZ (XYZ_r, XYZ_s, K_number, MakeName, DeviceName, DUTName, fnout)
            % Generate one-pager from XYZ
            
            XYZ_Reference = CPR.DeviceData;
            XYZ_Reference.setXYZ(XYZ_r);
            
            XYZ_Subject = CPR.DeviceData;
            XYZ_Subject.setXYZ(XYZ_s);
            
            cpr = CPR.ColorPerformanceReview(XYZ_Reference,XYZ_Subject);
            cpr.K_number = K_number;
            cpr.MakeName = MakeName;
            cpr.DeviceName = DeviceName;
            cpr.DUTName = DUTName;
            
            cpr.one_pager
            
            movefile('one_pager.png',fnout)
            
        end

    end

    %% Functions used in CPR Live Script (READEME_RST.mlx)
    %    function names started with "check_" 
    methods

        function check_visual_XYZ (obj)
            % Evaluate visual XYZ
            % Figure (a1)
            %   If the reproduced data were provided in the CIEXYZ color space (i.e., tristimulus),
            %   an optional Figure 8(a1) shows the color target by using CIE D65 as the reference white.
            %   The purpose is to check any excessive color shift caused by the endoscope light source and/or the device.
            
            [im0, im1, im2] = obj.evaluate_visual_XYZ_im;
            
            figure
            h_visual = gcf;
            
            subplot(1,3,1)
            imshow(im0)
            title('Original ColorChecker')
            
            subplot(1,3,2)
            imshow(im1)
            title('Reference Data, XYZ')
            
            subplot(1,3,3)
            imshow(im2)
            title('Reproduced Data, XYZ')
            
            set(h_visual,'Units','normalized','Position',[0 0 0.5 0.25]);
            obj.save_figure([48 267 676 229],'xyzx2.png')

            return
        end

        function check_Lab_visual (obj)
            % Visualize CIELAB data
            % The following charts show the simulated visual results when using the provided reference white. Use these charts to assess how the device would reproduce the ColorChecker.
        
            figure
            h_visual = gcf;
        
            ha0 = subplot(1,3,1);
            ha1 = subplot(1,3,2);
            ha2 = subplot(1,3,3);
            obj.evaluate_visual(ha0,ha1,ha2);
            set(h_visual,'Units','normalized','Position',[0 0 0.5 0.25]);

            obj.save_figure([48 267 676 229],'labx2.png')
        end
            
        function check_endoscopic_scene (obj,simulate_filename)
            % Visualize endoscopic scene
        
            figure
            im1 = imread(simulate_filename);
            ha1 = subplot(1,2,1);
            axes(ha1);
            imshow(im1);
            title('Ground Truth')
        
            im2 = obj.simulate_im (im1);
            ha2 = subplot(1,2,2);
            axes(ha2);
            imshow(im2);
            title('Reproduced')
            
            obj.save_figure([48 267 676 229],'samplex2.png')
        end
        
        function check_order (obj)
            % Order in lightness, chroma, and hue - 1D view
        
            figure
            hg = gcf;
            set(hg,'Visible','on')
        
            td = tiledlayout(4,1,'TileSpacing','tight');
            ha1 = nexttile(1); obj.order_1D_l6(ha1); title(sprintf('Lightness order of gray patches,{\\tau}=%.2f',obj.kendall_l6))
            ha2 = nexttile(2); obj.order_1D_l(ha2); title(sprintf('Lightness order of chromatic patches,{\\tau}=%.2f',obj.kendall_l))
            ha3 = nexttile(3); obj.order_1D_h(ha3); title(sprintf('Hue order of chromatic patches,{\\tau}=%.2f',obj.kendall_h))
            ha4 = nexttile(4); obj.order_1D_c(ha4); title(sprintf('Chroma order of chromatic patches,{\\tau}=%.2f',obj.kendall_c))
        
            obj.save_figure([31 380 673 478],'lch_1d.png')
        end
                
        function check_linearity (obj)
            % Order in lightness, chroma, and hue - 2D view
        
            figure
            h_order = gcf();
        
            ha1 = subplot(2,2,1);
            obj.evaluate_linearity_lightness_6(ha1); title('Lightness order of gray patches')
        
            ha2 = subplot(2,2,2);
            obj.evaluate_linearity_lightness(ha2); title('Lightness order of chromatic patches')
        
            ha3 = subplot(2,2,3);
            obj.evaluate_linearity_hue(ha3); title('Hue order of chromatic patches')
        
            ha4 = subplot(2,2,4);
            obj.evaluate_linearity_chroma(ha4); title('Chroma order of chromatic patches')
        
            obj.save_figure([436 321 834 777],'lch_2d.png')
        end

        function check_dE_from_truth (obj)
            % Absolute color errors in comparison with the ground truth
            % The left chart shows the per-patch color difference between the subject device and the ground truth. The right chart shows the bloxplot. Statistics (mean, std, min, median, and max) are provided in the titles.
            
            figure
            h_dE = gcf;
            ha1 = subplot(1,2,1);
            ha2 = subplot(1,2,2);
            obj.evaluate_dE_from_truth(ha1,ha2);
        
            %set(h_dE,'Units','normalized','Position',[0 0 0.5 0.25]);            
            obj.save_figure([48 267 676 229],'bar_de.png')
        end
          
        function check_color_transfer (obj)
            % Three-dimensional color transfer
        
            figure
            hg = gcf;
            set(hg,'Visible','on')
        
            td = tiledlayout(3,2,'TileSpacing','tight');
        
            ha1 = nexttile(td,1); obj.show3dquiver_select(obj.chromatic_patches,ha1); view([12 32]); axis equal; title('Chromatic patches')
            ha2 = nexttile(td,2); obj.show2dquiver_polar_select(obj.chromatic_patches,ha2); view([0 90]); axis([-100 100 -100 100 0 100]); title('Transfer in hue and chroma')
            ha3 = nexttile(td,3); obj.show3dquiver_select(obj.chromatic_patches,ha3); view([0 0]); axis square; axis([-100 100 -100 100 0 100]); title('Transfer in lightness and a*')
            ha4 = nexttile(td,4); obj.show3dquiver_select(obj.chromatic_patches,ha4); view([90 0]); axis square; axis([-100 100 -100 100 0 100]); title('Transfer in lightness and b*')
            ha5 = nexttile(td,5); obj.show3dquiver_select(obj.grayscale_patches,ha5); view([0 0]); axis square; axis([-100 100 -100 100 0 100]); title('Lightness of gray patches')
            ha6 = nexttile(td,6); obj.show2dquiver_polar_select(obj.grayscale_patches,ha6); view([0 90]); axis square; axis([-30 30 -30 30 0 100]); title('Chromaticity of gray patches')
        
            % capture the figure
            obj.save_figure([327 131 958 1075],'labx6.png')
        end

        function check_cce (obj)
            % Preservation of Color Contrast between Patches
        
            figure
            h_cce = gcf;
            ha1 = subplot(1,3,1); obj.evaluate_CCE(1,ha1);
            ha2 = subplot(1,3,2); obj.evaluate_CCE(2,ha2);
            ha3 = subplot(1,3,3); obj.evaluate_CCE(3,ha3);
            set(h_cce,'Units','normalized','Position',[0 0 1 1]);
        
            obj.save_figure([11 529 1227 556],'ccex3.png')
        end
        
    end
    
end


