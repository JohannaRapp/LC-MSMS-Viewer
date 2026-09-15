classdef LC_MSMS_Viewer_V2 < matlab.apps.AppBase
%LC_MSMS_VIEWER_V2  Browse targeted LC-MS/MS identifications.
%
%   Start the app by adding this folder and its subfolders to the MATLAB
%   path and running
%
%       LC_MSMS_Viewer_V2
%
%   The menus are
%
%       Load > Raw Data       run the whole pipeline on a set of mzXML files
%       Load > Analysed Data  reopen a data set saved earlier
%       Export > To Workspace copy the open data set to the base workspace
%
%   The table lists every targeted metabolite with its fold change and the
%   two quality-control verdicts. Selecting a row draws its chromatogram
%   against the median of the run, the MS1 spectrum inside the quadrupole
%   isolation window, and the three MS2 spectra with a star over every
%   fragment the CFM-ID prediction accounts for.
%
%   WHERE THE DATA LIVES
%
%   In the app's LcData property, not in a global variable. Three ways to
%   reach it:
%
%     Export > To Workspace   puts a copy in the base workspace as
%                             LC_MSMS_data, the name version 1 used
%
%     app = LC_MSMS_Viewer_V2;   start the app WITH an output argument,
%     d = app.LcData;            then read the property
%
%     load(savedFile)         a saved file holds one variable, also called
%                             LC_MSMS_data
%
%   Started without an output argument the handle is discarded, so use the
%   menu or the "app =" form to reach the data from the command line.
%
%   The whole pipeline also runs without the app:
%
%       lc = run_lc_msms_pipeline(rawFolder, hitList, predictions, outBase);
%
%   See also RUN_LC_MSMS_PIPELINE, LC_DEFAULT_PARAMS, PLOT_SAMPLE.
%
%   Authors: Hannes Link, Johanna Rapp.

    % --------------------------------------------------------------------
    % UI components
    % --------------------------------------------------------------------
    properties (Access = public)
        UIFigure             matlab.ui.Figure
        LoadMenu             matlab.ui.container.Menu
        LoadRawDataMenu      matlab.ui.container.Menu
        LoadAnalysedMenu     matlab.ui.container.Menu
        ExportMenu           matlab.ui.container.Menu
        ExportWorkspaceMenu  matlab.ui.container.Menu
        PolaritySwitch       matlab.ui.control.Switch
        PolaritySwitchLabel  matlab.ui.control.Label
        SampleTable          matlab.ui.control.Table
        ChromatogramAxes     matlab.ui.control.UIAxes
        IsolationAxes        matlab.ui.control.UIAxes
        Ms2Axes              matlab.ui.control.UIAxes
        StatusLabel          matlab.ui.control.Label
    end

    % --------------------------------------------------------------------
    % Application state
    % --------------------------------------------------------------------
    % Readable from outside so that "app = LC_MSMS_Viewer_V2; app.LcData"
    % works at the command line. Only the app writes to it. Version 1 used a
    % global, which made a second instance overwrite the first one's data.
    properties (GetAccess = public, SetAccess = private)
        LcData = []
    end

    properties (Access = private)
        SavedFilePath = ''
        Params
    end

    % --------------------------------------------------------------------
    % Menu callbacks
    % --------------------------------------------------------------------
    methods (Access = private)

        function LoadRawDataMenuSelected(app, ~)
            % Run the pipeline, then show the result.

            [lcData, savedFile] = run_lc_msms_pipeline('', '', '', '', app.Params);

            if isempty(lcData)
                return      % a dialog was cancelled
            end

            app.LcData = lcData;
            app.SavedFilePath = savedFile;
            app.refreshAll();
        end

        function LoadAnalysedMenuSelected(app, ~)
            % Reopen a data set saved by an earlier run.

            [fileName, folderPath] = uigetfile('*.mat', 'Select an analysed LC-MS/MS data set');
            if isequal(fileName, 0)
                return
            end

            fullName = fullfile(folderPath, fileName);
            loaded = load(fullName, 'LC_MSMS_data');

            if ~isfield(loaded, 'LC_MSMS_data')
                uialert(app.UIFigure, ...
                    sprintf('%s holds no variable called LC_MSMS_data.', fileName), ...
                    'Not an LC-MS/MS data set');
                return
            end

            app.LcData = loaded.LC_MSMS_data;
            app.SavedFilePath = fullName;

            % Work with the settings the data were produced with.
            if isfield(app.LcData, 'params') && ~isempty(app.LcData(1).params)
                app.Params = app.LcData(1).params;
            end

            app.refreshAll();
        end

        function ExportWorkspaceMenuSelected(app, ~)
            if ~app.hasData()
                uialert(app.UIFigure, 'Load a data set first.', 'No data');
                return
            end

            assignin('base', 'LC_MSMS_data', app.LcData);

            fprintf(['LC_MSMS_data is now in your base workspace. It is a ' ...
                'copy: editing it there does not change what the app holds.\n']);

            uialert(app.UIFigure, ...
                ['The variable LC_MSMS_data is now in your base workspace. ' ...
                 'It is a copy, so editing it there does not change what ' ...
                 'the app holds.'], 'Sent to workspace', 'Icon', 'success');
        end

    end

    % --------------------------------------------------------------------
    % Table and control callbacks
    % --------------------------------------------------------------------
    methods (Access = private)

        function SampleTableCellSelection(app, event)
            if isempty(event.Indices) || ~app.hasData()
                return
            end

            modeData = app.LcData(app.currentMode());
            row = event.Indices(1);

            if row > height(modeData.samples)
                return
            end

            plot_sample(modeData, row, app.ChromatogramAxes, ...
                app.IsolationAxes, app.Ms2Axes);
        end

        function PolaritySwitchValueChanged(app, ~)
            app.clearPlots();
            app.refreshTable();
        end

    end

    % --------------------------------------------------------------------
    % Display helpers
    % --------------------------------------------------------------------
    methods (Access = private)

        function tf = hasData(app)
            tf = ~isempty(app.LcData) && isfield(app.LcData, 'samples') ...
                && ~isempty(app.LcData(app.currentMode()).samples);
        end

        function mode = currentMode(app)
            % 1 is positive polarity, 2 negative.
            if strcmp(app.PolaritySwitch.Value, 'pos')
                mode = 1;
            else
                mode = 2;
            end
        end

        function refreshAll(app)
            app.clearPlots();
            app.refreshTable();

            if isempty(app.SavedFilePath)
                app.StatusLabel.Text = '(not saved)';
            else
                app.StatusLabel.Text = app.SavedFilePath;
            end
        end

        function refreshTable(app)
            if ~app.hasData()
                app.SampleTable.Data = table();
                return
            end

            modeData = app.LcData(app.currentMode());

            app.SampleTable.Data = table( ...
                modeData.samples.abbreviation, ...
                modeData.samples.gene, ...
                modeData.foldChange, ...
                [modeData.purity.intensityVerdict].', ...
                [modeData.purity.sidePeakVerdict].', ...
                'VariableNames', {'Metabolite', 'Gene', 'FoldChange', ...
                                  'QC_intensity', 'QC_sidePeaks'});
        end

        function clearPlots(app)
            cla(app.ChromatogramAxes);
            cla(app.IsolationAxes);
            for k = 1:numel(app.Ms2Axes)
                cla(app.Ms2Axes(k));
            end
        end

    end

    % --------------------------------------------------------------------
    % Component initialisation
    % --------------------------------------------------------------------
    methods (Access = private)

        function createComponents(app)

            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [1 1 1];
            app.UIFigure.Position = [100 100 1000 620];
            app.UIFigure.Name = 'LC-MS/MS Viewer V2';

            % --- menus ---
            app.LoadMenu = uimenu(app.UIFigure);
            app.LoadMenu.Text = 'Load';

            app.LoadRawDataMenu = uimenu(app.LoadMenu);
            app.LoadRawDataMenu.Text = 'Raw Data';
            app.LoadRawDataMenu.MenuSelectedFcn = ...
                createCallbackFcn(app, @LoadRawDataMenuSelected, true);

            app.LoadAnalysedMenu = uimenu(app.LoadMenu);
            app.LoadAnalysedMenu.Text = 'Analysed Data';
            app.LoadAnalysedMenu.MenuSelectedFcn = ...
                createCallbackFcn(app, @LoadAnalysedMenuSelected, true);

            app.ExportMenu = uimenu(app.UIFigure);
            app.ExportMenu.Text = 'Export';

            app.ExportWorkspaceMenu = uimenu(app.ExportMenu);
            app.ExportWorkspaceMenu.Text = 'To Workspace (LC_MSMS_data)';
            app.ExportWorkspaceMenu.MenuSelectedFcn = ...
                createCallbackFcn(app, @ExportWorkspaceMenuSelected, true);

            % --- polarity ---
            app.PolaritySwitchLabel = uilabel(app.UIFigure);
            app.PolaritySwitchLabel.Position = [24 578 50 22];
            app.PolaritySwitchLabel.Text = 'Polarity';

            app.PolaritySwitch = uiswitch(app.UIFigure, 'slider');
            app.PolaritySwitch.Items = {'pos', 'neg'};
            app.PolaritySwitch.Value = 'neg';
            app.PolaritySwitch.Position = [92 580 45 20];
            app.PolaritySwitch.ValueChangedFcn = ...
                createCallbackFcn(app, @PolaritySwitchValueChanged, true);

            % --- table ---
            app.SampleTable = uitable(app.UIFigure);
            app.SampleTable.RowName = {};
            app.SampleTable.Position = [24 300 430 260];
            app.SampleTable.CellSelectionCallback = ...
                createCallbackFcn(app, @SampleTableCellSelection, true);

            % --- plots ---
            app.ChromatogramAxes = uiaxes(app.UIFigure);
            app.ChromatogramAxes.Position = [480 320 490 250];

            app.IsolationAxes = uiaxes(app.UIFigure);
            app.IsolationAxes.Position = [24 30 430 250];

            for k = 1:3
                axesHandle = uiaxes(app.UIFigure);
                axesHandle.Position = [480 + (k-1)*170, 30, 165, 250];
                app.Ms2Axes(k) = axesHandle;
            end

            % --- status ---
            app.StatusLabel = uilabel(app.UIFigure);
            app.StatusLabel.FontSize = 9;
            app.StatusLabel.Position = [24 566 946 16];
            app.StatusLabel.Text = '';

            app.UIFigure.Visible = 'on';
        end

    end

    % --------------------------------------------------------------------
    % Construction and destruction
    % --------------------------------------------------------------------
    methods (Access = public)

        function app = LC_MSMS_Viewer_V2()

            % Make the pipeline reachable even when only this folder is on
            % the MATLAB path.
            appFolder = fileparts(mfilename('fullpath'));
            addpath(appFolder, fullfile(appFolder, 'Functions'));

            app.Params = lc_default_params();

            createComponents(app);
            registerApp(app, app.UIFigure);

            if nargout == 0
                clear app
            end
        end

        function delete(app)
            delete(app.UIFigure);
        end

    end

end
