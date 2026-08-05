classdef AP_Analysis_APP < matlab.apps.AppBase

    properties (Access = public)
        UIFigure
        MainGrid
        HeaderPanel
        HeaderGrid
        TitleLabel
        StatusLamp
        StatusLabel
        ModeLabel
        ModeDropDown
        ScopeLabel
        ScopeDropDown
        InputLabel
        InputPathEdit
        BrowseInputButton
        ReuseSourceLabel
        ReuseSourceEdit
        BrowseReuseSourceButton
        DetectReuseSourceButton
        PresetLabel
        PresetDropDown
        RecordAnalysisPanel
        RecordAnalysisGrid
        RecordAverageCheckBox
        RecordAverageOnlyCheckBox
        LoadRequestButton
        SaveRequestButton
        ValidateButton
        RunButton
        OpenOutputButton
        MainTabs
        WorkflowTab
        WorkflowGrid
        WorkflowHelpLabel
        SectionTable
        BasicTab
        BasicGrid
        BasicSwitchLabel
        BasicSwitchTable
        BasicValueLabel
        BasicValueTable
        AdvancedTab
        AdvancedGrid
        AdvancedSwitchLabel
        AdvancedSwitchTable
        AdvancedValueLabel
        AdvancedValueTable
        ExpertTab
        ExpertGrid
        ExpertSwitchLabel
        ExpertSwitchTable
        ExpertValueLabel
        ExpertValueTable
        LogDivider
        LogPanel
        LogGrid
        LogHelpLabel
        LogTextArea
    end

    properties (Access = private)
        ProjectRoot string = ""
        Request struct = struct()
        ParameterCatalog struct = struct([])
        SectionCatalog struct = struct([])
        CurrentOutputPath string = ""
        LogLines string = strings(0, 1)
        IsBusy logical = false
        SessionLogger = []
        LogResizeActive logical = false
        LogResizeStartY double = NaN
        LogResizeStartHeight double = 320
        LogResizeAvailableHeight double = NaN
        AuthoritativePreflight struct = struct()
    end

    methods (Access = private)

        function startupFcn(app)
            app.ProjectRoot = string(fileparts(which('AP_Analysis_APP')));
            if strlength(app.ProjectRoot) == 0
                app.ProjectRoot = string(pwd);
            end
            startup_file = fullfile(app.ProjectRoot, 'AAA_startup.m');
            if isfile(startup_file)
                addpath(app.ProjectRoot, '-begin');
                AAA_startup();
            else
                error('AAA:App:StartupMissing', ...
                    'AAA_startup.m was not found beside AP_Analysis_APP.mlapp.');
            end

            session_log_dir = fullfile(app.ProjectRoot, 'runtime', 'logs');
            if ~isfolder(session_log_dir), mkdir(session_log_dir); end
            session_log_file = fullfile(session_log_dir, ...
                char("AAA_app_session_" + string(datetime('now', ...
                'Format', 'yyyy-MM-dd_HH-mm-ss-SSS')) + ".log"));
            app.SessionLogger = aaa.io.AnalysisLogger(session_log_file, false);
            app.appendLog("AAA 已启动。所有运算由 aaa.run_analysis3 分发。");
            app.appendLog("App session log: " + string(session_log_file));
            app.resetRequest("dual", "single");
        end

        function resetRequest(app, mode, scope)
            app.Request = aaa.schema.default_request(mode, scope);
            app.CurrentOutputPath = "";
            app.AuthoritativePreflight = struct();
            app.syncRequestToUI();
            app.appendLog("已载入 " + upper(string(mode)) + ...
                " / " + string(scope) + " 默认 request；模式切换会重置参数。 ");
        end

        function syncRequestToUI(app)
            [mode, found] = app.getPath(app.Request, 'mode');
            if ~found, mode = "dual"; end
            [scope, found] = app.getPath(app.Request, 'scope');
            if ~found, scope = "single"; end
            mode = lower(string(mode));
            scope = lower(string(scope));

            app.ModeDropDown.Value = mode;
            app.ScopeDropDown.Value = scope;

            [input_path, found] = app.getPath(app.Request, 'input_path');
            if ~found, input_path = ""; end
            app.InputPathEdit.Value = char(string(input_path));

            [reuse_path, found] = app.getPath(app.Request, 'workflow.source_results_path');
            if ~found, reuse_path = ""; end
            app.ReuseSourceEdit.Value = char(string(reuse_path));

            [preset, found] = app.getPath(app.Request, 'workflow.preset');
            if ~found || strlength(string(preset)) == 0, preset = "full"; end
            preset = lower(string(preset));
            if ~ismember(preset, string(app.PresetDropDown.ItemsData))
                preset = "custom";
            end
            app.PresetDropDown.Value = preset;

            [record_average, found] = app.getPath(app.Request, 'record.record_average');
            if ~found, record_average = false; end
            [record_average_only, found] = app.getPath(app.Request, 'record.record_average_only');
            if ~found, record_average_only = false; end
            is_record = scope == "record";
            app.RecordAnalysisPanel.Visible = app.onOff(is_record);
            app.RecordAverageCheckBox.Value = logical(record_average);
            app.RecordAverageOnlyCheckBox.Value = logical(record_average_only);
            app.RecordAverageOnlyCheckBox.Enable = app.onOff(is_record && logical(record_average));

            app.ParameterCatalog = aaa.schema.parameter_catalog(mode, scope);
            app.SectionCatalog = aaa.schema.section_catalog(mode);
            app.ensureResolvedSections();
            app.refreshSectionTable();
            app.refreshParameterTables();
            app.StatusLabel.Text = sprintf('%s / %s · request ready', ...
                upper(char(mode)), char(scope));
        end

        function ensureResolvedSections(app)
            [preset, found] = app.getPath(app.Request, 'workflow.preset');
            if ~found || strlength(string(preset)) == 0, preset = "full"; end
            [sections, found] = app.getPath(app.Request, 'workflow.sections');
            if ~found || ~isstruct(sections)
                sections = struct();
            end
            if string(preset) == "custom"
                names = [app.SectionCatalog.name];
                for idx = 1:numel(names)
                    name = char(names(idx));
                    if ~isfield(sections, name)
                        sections.(name) = app.SectionCatalog(idx).allowed_actions(1);
                    end
                end
            end
            [~, ~, ~, resolved] = aaa.schema.resolve_sections( ...
                app.Request.mode, preset, sections);
            app.Request = app.setPath(app.Request, 'workflow.sections', resolved);
        end

        function refreshSectionTable(app)
            catalog = app.SectionCatalog;
            n = numel(catalog);
            names = strings(n, 1);
            labels = strings(n, 1);
            actions = strings(n, 1);
            statuses = repmat("pending", n, 1);
            descriptions = strings(n, 1);
            [sections, ~] = app.getPath(app.Request, 'workflow.sections');
            old_data = app.SectionTable.Data;
            for idx = 1:n
                names(idx) = catalog(idx).name;
                labels(idx) = catalog(idx).ui_label;
                descriptions(idx) = catalog(idx).description;
                name = char(names(idx));
                if isstruct(sections) && isfield(sections, name)
                    actions(idx) = string(sections.(name));
                else
                    actions(idx) = catalog(idx).allowed_actions(1);
                end
                if istable(old_data) && ismember('Section', old_data.Properties.VariableNames)
                    old_idx = find(string(old_data.Section) == names(idx), 1);
                    if ~isempty(old_idx) && ismember('Status', old_data.Properties.VariableNames)
                        statuses(idx) = string(old_data.Status(old_idx));
                    end
                end
            end
            actions = categorical(actions, ["run", "reuse", "skip"]);
            app.SectionTable.Data = table(names, labels, actions, statuses, descriptions, ...
                'VariableNames', {'Section','Label','Action','Status','Description'});
        end

        function refreshParameterTables(app)
            app.populateParameterLevel("Basic");
            app.populateParameterLevel("Advanced");
            app.populateParameterLevel("Expert");
        end

        function populateParameterLevel(app, level)
            catalog = app.ParameterCatalog;
            if isempty(catalog)
                logical_catalog = catalog;
                value_catalog = catalog;
            else
                keep = [catalog.ui_level] == level & ...
                    ~ismember([catalog.key], ["input_path", "workflow.preset", ...
                    "workflow.source_results_path","record.record_average", ...
                    "record.record_average_only"]);
                level_catalog = catalog(keep);
                logical_catalog = level_catalog([level_catalog.value_type] == "logical");
                value_catalog = level_catalog([level_catalog.value_type] ~= "logical");
            end

            switch_data = app.makeSwitchTableData(logical_catalog);
            value_data = app.makeValueTableData(value_catalog);
            switch level
                case "Basic"
                    app.BasicSwitchTable.Data = switch_data;
                    app.BasicValueTable.Data = value_data;
                case "Advanced"
                    app.AdvancedSwitchTable.Data = switch_data;
                    app.AdvancedValueTable.Data = value_data;
                case "Expert"
                    app.ExpertSwitchTable.Data = switch_data;
                    app.ExpertValueTable.Data = value_data;
            end
        end

        function data = makeSwitchTableData(app, catalog)
            n = numel(catalog);
            parameter = strings(n, 1);
            enabled = false(n, 1);
            group = strings(n, 1);
            section = strings(n, 1);
            key = strings(n, 1);
            description = strings(n, 1);
            for idx = 1:n
                entry = catalog(idx);
                parameter(idx) = app.parameterDisplayLabel(entry);
                key(idx) = entry.key;
                group(idx) = entry.group;
                section(idx) = entry.section;
                description(idx) = app.decorateDescription(entry);
                [value, found] = app.getPath(app.Request, entry.key);
                if ~found, value = entry.default; end
                enabled(idx) = logical(value);
            end
            data = table(parameter, enabled, group, section, key, description, ...
                'VariableNames', {'Parameter','Enabled','Group','Section','Key','Description'});
        end

        function data = makeValueTableData(app, catalog)
            n = numel(catalog);
            parameter = strings(n, 1);
            value = strings(n, 1);
            units = strings(n, 1);
            allowed = strings(n, 1);
            group = strings(n, 1);
            section = strings(n, 1);
            key = strings(n, 1);
            description = strings(n, 1);
            for idx = 1:n
                entry = catalog(idx);
                parameter(idx) = app.parameterDisplayLabel(entry);
                key(idx) = entry.key;
                group(idx) = entry.group;
                section(idx) = entry.section;
                units(idx) = entry.units;
                if isempty(entry.choices)
                    allowed(idx) = "";
                else
                    allowed(idx) = strjoin(entry.choices, " | ");
                end
                description(idx) = app.decorateDescription(entry);
                [request_value, found] = app.getPath(app.Request, entry.key);
                if ~found, request_value = entry.default; end
                value(idx) = app.valueToText(request_value);
            end
            data = table(parameter, value, units, allowed, group, section, key, description, ...
                'VariableNames', {'Parameter','Value','Units','Allowed','Group','Section','Key','Description'});
        end

        function text_value = decorateDescription(~, entry)
            text_value = string(entry.description);
            if strlength(entry.visible_when) > 0
                text_value = text_value + " [相关条件: " + entry.visible_when + "]";
            end
            if strlength(entry.legacy_source) > 0
                text_value = text_value + " [来源: " + entry.legacy_source + "]";
            end
            if isfield(entry, 'default_rule') && strlength(entry.default_rule) > 0
                text_value = text_value + " [默认规则: " + entry.default_rule + "]";
            end
            if isfield(entry, 'editable') && ~entry.editable
                text_value = "[只读] " + text_value;
            end
        end

        function label = parameterDisplayLabel(~, entry)
            label = string(entry.ui_label);
            if isfield(entry, 'editable') && ~entry.editable
                label = "[只读] " + label;
            end
        end

        function ModeDropDownValueChanged(app, ~)
            if app.IsBusy, return; end
            app.appendLog("[CONTROL] Mode changed: "+string(app.Request.mode)+ ...
                " -> "+string(app.ModeDropDown.Value)+".");
            app.resetRequest(string(app.ModeDropDown.Value), ...
                string(app.ScopeDropDown.Value));
            app.appendLog("[PREFLIGHT] Previous validation cleared by Mode change.");
            drawnow limitrate;
        end

        function ScopeDropDownValueChanged(app, ~)
            if app.IsBusy, return; end
            app.appendLog("[CONTROL] Scope changed: "+string(app.Request.scope)+ ...
                " -> "+string(app.ScopeDropDown.Value)+".");
            app.resetRequest(string(app.ModeDropDown.Value), ...
                string(app.ScopeDropDown.Value));
            app.appendLog("[PREFLIGHT] Previous validation cleared by Scope change.");
            drawnow limitrate;
        end

        function InputPathEditValueChanged(app, ~)
            old_path=string(app.Request.input_path);
            new_path=string(app.InputPathEdit.Value);
            app.Request = app.setPath(app.Request, 'input_path', ...
                new_path);
            app.appendLog("[CONTROL] Input path changed: "+old_path+" -> "+new_path+".");
            app.invalidatePreflight("Input path change");
            app.autoDetectScopeFromInput();
            app.quickInspectAndReport();
        end

        function BrowseInputButtonPushed(app, ~)
            app.logButtonClick("Browse input");
            initial_path = app.InputPathEdit.Value;
            if ~isfolder(initial_path)
                initial_path = char(app.ProjectRoot);
            end
            selected = uigetdir(initial_path, '选择 Cycle 或 Record 输入目录');
            if isequal(selected, 0)
                app.appendLog("[UI] Browse input cancelled.");
                return;
            end
            app.InputPathEdit.Value = selected;
            app.InputPathEditValueChanged([]);
        end

        function ReuseSourceEditValueChanged(app, ~)
            old_path=string(app.Request.workflow.source_results_path);
            new_path=string(app.ReuseSourceEdit.Value);
            app.Request = app.setPath(app.Request, ...
                'workflow.source_results_path',new_path);
            app.appendLog("[CONTROL] Manual Reuse source changed: "+old_path+" -> "+new_path+".");
            app.invalidatePreflight("Reuse source change");
            if strlength(strtrim(new_path))==0
                app.appendLog("[QUICK CHECK] Reuse source is blank; automatic detection remains pending.");
            elseif isfile(new_path)||isfolder(new_path)
                app.appendLog("[QUICK CHECK] Reuse source path exists; compatibility remains pending.");
            else
                app.appendLog("[WARN] Reuse source path does not exist: "+new_path);
            end
            drawnow limitrate;
        end

        function BrowseReuseSourceButtonPushed(app, ~)
            app.logButtonClick("Browse Reuse source");
            initial_path = app.ReuseSourceEdit.Value;
            if ~isfolder(initial_path), initial_path = app.InputPathEdit.Value; end
            if ~isfolder(initial_path), initial_path = char(app.ProjectRoot); end
            selected = uigetdir(initial_path, 'Select reuse result or search root');
            if isequal(selected,0)
                app.appendLog("[UI] Browse Reuse source cancelled.");
                return;
            end
            app.ReuseSourceEdit.Value = selected;
            app.ReuseSourceEditValueChanged([]);
        end

        function DetectReuseSourceButtonPushed(app, ~)
            app.logButtonClick("Detect Reuse source");
            app.syncTopControlsToRequest();
            app.performAuthoritativePreflight("Detect reuse");
        end

        function PresetDropDownValueChanged(app, ~)
            if app.IsBusy, return; end
            old_preset=string(app.Request.workflow.preset);
            preset = string(app.PresetDropDown.Value);
            app.Request = app.setPath(app.Request, 'workflow.preset', preset);
            if preset == "custom"
                [current, found] = app.getPath(app.Request, 'workflow.sections');
                if ~found, current = struct(); end
                [~, ~, ~, resolved] = aaa.schema.resolve_sections( ...
                    app.Request.mode, preset, current);
            else
                [~, ~, ~, resolved] = aaa.schema.resolve_sections( ...
                    app.Request.mode, preset, struct());
            end
            app.Request = app.setPath(app.Request, 'workflow.sections', resolved);
            app.refreshSectionTable();
            app.appendLog("[CONTROL] Preset changed: "+old_preset+" -> "+preset+".");
            app.invalidatePreflight("Preset change");
            drawnow limitrate;
        end

        function RecordAverageCheckBoxValueChanged(app, ~)
            if app.IsBusy, return; end
            old_value = logical(app.Request.record.record_average);
            new_value = logical(app.RecordAverageCheckBox.Value);
            app.Request.record.record_average = new_value;
            if ~new_value
                app.Request.record.record_average_only = false;
                app.RecordAverageOnlyCheckBox.Value = false;
                app.RecordAverageOnlyCheckBox.Enable = 'off';
            else
                app.RecordAverageOnlyCheckBox.Enable = 'on';
            end
            app.appendLog("Record analysis changed: " + app.onOff(old_value) + ...
                " -> " + app.onOff(new_value) + ".");
            app.invalidatePreflight("Record analysis switch change");
            drawnow limitrate;
        end

        function RecordAverageOnlyCheckBoxValueChanged(app, ~)
            if app.IsBusy, return; end
            old_value = logical(app.Request.record.record_average_only);
            new_value = logical(app.RecordAverageOnlyCheckBox.Value);
            if new_value && ~app.RecordAverageCheckBox.Value
                app.RecordAverageCheckBox.Value = true;
                app.Request.record.record_average = true;
            end
            app.Request.record.record_average_only = new_value;
            app.appendLog("Record source changed: existing Cycle results only " + ...
                app.onOff(old_value) + " -> " + app.onOff(new_value) + ".");
            app.invalidatePreflight("Record-average source switch change");
            drawnow limitrate;
        end

        function SectionTableCellEdit(app, event)
            if app.IsBusy || isempty(event.Indices) || event.Indices(2) ~= 3
                return;
            end
            row = event.Indices(1);
            data = app.SectionTable.Data;
            action = lower(string(event.NewData));
            if ~ismember(action, ["run", "reuse", "skip"])
                app.refreshSectionTable();
                app.appendLog("[ERROR] Invalid Section action entered: "+action+".");
                drawnow limitrate;
                uialert(app.UIFigure, 'Section action 必须是 run、reuse 或 skip。', ...
                    '无效的 section action');
                return;
            end
            name = char(string(data.Section(row)));
            catalog_idx = find([app.SectionCatalog.name] == string(name), 1);
            if isempty(catalog_idx)
                app.refreshSectionTable();
                uialert(app.UIFigure, sprintf('未知 section: %s。', name), ...
                    '无效的 section');
                return;
            end
            allowed_actions = app.SectionCatalog(catalog_idx).allowed_actions;
            if ~ismember(action, allowed_actions)
                app.refreshSectionTable();
                app.appendLog("[ERROR] Section "+string(name)+" does not allow action "+action+".");
                drawnow limitrate;
                uialert(app.UIFigure, ...
                    sprintf('%s 允许的 action: %s。', name, ...
                    strjoin(allowed_actions, ', ')), ...
                    '该 section 不支持此 action');
                return;
            end
            [sections, ~] = app.getPath(app.Request, 'workflow.sections');
            old_action=string(sections.(name));
            old_preset=string(app.Request.workflow.preset);
            sections.(name) = action;
            app.Request = app.setPath(app.Request, 'workflow.sections', sections);
            app.Request = app.setPath(app.Request, 'workflow.preset', "custom");
            app.PresetDropDown.Value = "custom";
            app.appendLog("[CONTROL] Section "+string(name)+" action changed: "+ ...
                old_action+" -> "+action+".");
            if old_preset~="custom"
                app.appendLog("[CONTROL] Preset changed: "+old_preset+" -> custom.");
            end
            app.invalidatePreflight("Section "+string(name)+" action change");
            drawnow limitrate;
        end

        function SwitchTableCellEdit(app, source, event)
            if app.IsBusy || isempty(event.Indices) || event.Indices(2) ~= 2
                return;
            end
            row = event.Indices(1);
            data = source.Data;
            key = string(data.Key(row));
            match = find([app.ParameterCatalog.key] == key, 1);
            if ~isempty(match) && isfield(app.ParameterCatalog(match), 'editable') ...
                    && ~app.ParameterCatalog(match).editable
                app.refreshParameterTables();
                uialert(app.UIFigure, '该值由 backend 自动决定，只读展示。', ...
                    '只读参数');
                return;
            end
            [old_value,found]=app.getPath(app.Request,key);
            if ~found,old_value=[];end
            new_value=logical(event.NewData);
            app.Request = app.setPath(app.Request,key,new_value);
            app.appendLog(key+" changed: "+app.valueToText(old_value)+ ...
                " -> "+app.valueToText(new_value)+".");
            app.invalidatePreflight(key+" change");
            drawnow limitrate;
        end

        function ValueTableCellEdit(app, source, event)
            if app.IsBusy || isempty(event.Indices) || event.Indices(2) ~= 2
                return;
            end
            row = event.Indices(1);
            data = source.Data;
            key = string(data.Key(row));
            match = find([app.ParameterCatalog.key] == key, 1);
            if isempty(match)
                app.refreshParameterTables();
                return;
            end
            entry = app.ParameterCatalog(match);
            if isfield(entry, 'editable') && ~entry.editable
                app.refreshParameterTables();
                uialert(app.UIFigure, ...
                    '该值由 backend 自动决定，只读展示，不能手工覆盖。', ...
                    '只读参数');
                return;
            end
            try
                value = app.parseTextValue(event.NewData, entry);
                [old_value,found]=app.getPath(app.Request,key);
                if ~found,old_value=[];end
                app.Request = app.setPath(app.Request, key, value);
                app.appendLog(key+" changed: "+app.valueToText(old_value)+ ...
                    " -> "+app.valueToText(value)+".");
                app.invalidatePreflight(key+" change");
                drawnow limitrate;
            catch err
                app.refreshParameterTables();
                uialert(app.UIFigure, err.message, '参数格式错误');
            end
        end

        function LoadRequestButtonPushed(app, ~)
            app.logButtonClick("Load request");
            [file, folder] = uigetfile('*.mat', '载入 AAA request');
            if isequal(file, 0)
                app.appendLog("[UI] Load request cancelled.");
                return;
            end
            request_file=string(fullfile(folder,file));
            try
                loaded=load(request_file);
            catch err
                app.appendLog("[ERROR] Failed to load request: "+string(err.message));
                uialert(app.UIFigure,err.message,'request load failed');
                return;
            end
            if ~isfield(loaded, 'request') || ~isstruct(loaded.request)
                app.appendLog("[ERROR] Selected MAT file has no struct named request: "+request_file);
                uialert(app.UIFigure, 'MAT 文件必须包含结构体变量 request。', ...
                    '无法载入 request');
                return;
            end
            try
                [normalized, report] = aaa.schema.validate_request(loaded.request);
                app.Request = normalized;
                app.syncRequestToUI();
                app.appendLog("Request loaded from: "+string(fullfile(folder,file))+".");
                app.invalidatePreflight("Loaded request");
                app.quickInspectAndReport();
                app.reportValidation(report, "载入 request");
            catch err
                app.appendLog("[ERROR] Loaded request is invalid: "+string(err.message));
                uialert(app.UIFigure, err.message, 'request 无效');
            end
        end

        function SaveRequestButtonPushed(app, ~)
            app.logButtonClick("Save request");
            app.syncTopControlsToRequest();
            [file, folder] = uiputfile('*.mat', '保存 AAA request', ...
                sprintf('AAA_%s_%s_request.mat', app.Request.mode, app.Request.scope));
            if isequal(file, 0)
                app.appendLog("[UI] Save request cancelled.");
                return;
            end
            request = app.Request; %#ok<NASGU>
            request_file=string(fullfile(folder,file));
            try
                save(request_file, 'request', '-v7.3');
                app.appendLog("Request saved: " + request_file);
            catch err
                app.appendLog("[ERROR] Failed to save request: "+string(err.message));
                uialert(app.UIFigure,err.message,'request save failed');
            end
        end

        function ValidateButtonPushed(app, ~)
            app.logButtonClick("Validate request: "+upper(string(app.ModeDropDown.Value))+ ...
                " / "+string(app.ScopeDropDown.Value));
            app.syncTopControlsToRequest();
            try
                preflight=app.performAuthoritativePreflight("Validate request");
                if ~preflight.valid
                    error('AAA:App:InvalidPreflight','%s', ...
                        app.preflightErrorText(preflight));
                end
                app.StatusLabel.Text='request valid';
            catch err
                app.StatusLabel.Text = 'validation failed';
                app.appendLog("[ERROR] " + string(err.message));
                uialert(app.UIFigure, err.message, '验证失败');
            end
        end

        function RunButtonPushed(app, ~)
            app.logButtonClick("Run AAA: "+upper(string(app.ModeDropDown.Value))+ ...
                " / "+string(app.ScopeDropDown.Value));
            if app.IsBusy
                app.appendLog("[UI] Run ignored because AAA is already busy.");
                return;
            end
            app.syncTopControlsToRequest();
            try
                [preflight,reused]=app.preflightForRun();
                if reused
                    app.appendLog("[PREFLIGHT] Run found a matching validation record; deep scan skipped.");
                end
                if ~preflight.valid
                    if reused
                        app.appendLog("[PREFLIGHT] The matching record is failed; Run remains blocked without rescanning.");
                    end
                    app.appendLog("[ERROR] Run blocked by authoritative preflight.");
                    uialert(app.UIFigure,app.preflightErrorText(preflight), ...
                        'Authoritative preflight failed');
                    return;
                end
                run_request=preflight.resolved_request;
            catch err
                app.appendLog("[ERROR] Run preflight failed: "+string(err.message));
                uialert(app.UIFigure, err.message, '运行前验证失败');
                return;
            end

            app.setBusy(true);
            cleanup = onCleanup(@() app.setBusy(false)); %#ok<NASGU>
            app.resetSectionStatuses();
            app.appendLog("开始运行 request。 ");
            runtime = struct();
            runtime.notify = @(event) app.onRuntimeEvent(event);
            runtime.event_callback = runtime.notify;
            runtime.app_root = app.ProjectRoot;
            runtime.preflight = preflight;
            if ~isempty(app.SessionLogger) && isvalid(app.SessionLogger)
                runtime.logger = app.SessionLogger;
            end
            try
                receipt = aaa.run_analysis3(run_request, runtime);
                app.captureReceipt(receipt);
                app.StatusLabel.Text = 'completed';
                app.appendLog("运行完成。 ");
            catch err
                app.StatusLabel.Text = 'failed';
                app.appendLog("[ERROR] " + string(getReport(err, 'extended', ...
                    'hyperlinks', 'off')));
                uialert(app.UIFigure, err.message, 'AAA 运行失败');
            end
        end

        function OpenOutputButtonPushed(app, ~)
            app.logButtonClick("Open output");
            output_path = app.CurrentOutputPath;
            if strlength(output_path) == 0
                [candidate, found] = app.getPath(app.Request, 'workflow.output_path');
                if found, output_path = string(candidate); end
            end
            if strlength(output_path) == 0 || ~isfolder(output_path)
                app.appendLog("[WARN] Open output failed: no existing output folder is available.");
                uialert(app.UIFigure, '尚无可打开的输出目录。', '输出目录');
                return;
            end
            app.appendLog("Opening output folder: "+output_path);
            winopen(char(output_path));
        end

        function syncTopControlsToRequest(app)
            app.Request = app.setPath(app.Request, 'mode', ...
                lower(string(app.ModeDropDown.Value)));
            app.Request = app.setPath(app.Request, 'scope', ...
                lower(string(app.ScopeDropDown.Value)));
            app.Request = app.setPath(app.Request, 'input_path', ...
                string(app.InputPathEdit.Value));
            app.Request = app.setPath(app.Request, 'workflow.source_results_path', ...
                string(app.ReuseSourceEdit.Value));
            app.Request = app.setPath(app.Request, 'workflow.preset', ...
                lower(string(app.PresetDropDown.Value)));
            if lower(string(app.ScopeDropDown.Value)) == "record"
                app.Request = app.setPath(app.Request, 'record.record_average', ...
                    logical(app.RecordAverageCheckBox.Value));
                app.Request = app.setPath(app.Request, 'record.record_average_only', ...
                    logical(app.RecordAverageOnlyCheckBox.Value));
            end
        end

        function reportValidation(app, report, context)
            lines = app.reportLines(report);
            if isempty(lines)
                lines = context + ": 验证完成。";
            end
            for idx = 1:numel(lines)
                app.appendLog(lines(idx));
            end
            if isfield(report, 'valid') && report.valid
                app.StatusLabel.Text = 'request valid';
            else
                app.StatusLabel.Text = 'request invalid';
            end
        end

        function lines = reportLines(~, report)
            lines = strings(0, 1);
            if isfield(report, 'errors') && ~isempty(report.errors)
                values = string(report.errors(:));
                lines = [lines; "[ERROR] " + values]; %#ok<AGROW>
            end
            if isfield(report, 'warnings') && ~isempty(report.warnings)
                values = string(report.warnings(:));
                lines = [lines; "[WARN] " + values]; %#ok<AGROW>
            end
            if isfield(report, 'valid') && report.valid && isempty(lines)
                lines = "[OK] request valid";
            end
        end

        function onRuntimeEvent(app, event)
            try
                if isstruct(event)
                    message = "";
                    if isfield(event, 'message'), message = string(event.message); end
                    if strlength(message) == 0 && isfield(event, 'type')
                        message = string(event.type);
                    end
                    prefix = strings(0, 1);
                    if isfield(event, 'cycle') && strlength(string(event.cycle)) > 0
                        prefix(end+1) = "cycle=" + string(event.cycle); %#ok<AGROW>
                    end
                    if isfield(event, 'section') && strlength(string(event.section)) > 0
                        prefix(end+1) = "section=" + string(event.section); %#ok<AGROW>
                        status = "running";
                        if isfield(event, 'state') && strlength(string(event.state)) > 0
                            status = string(event.state);
                        elseif isfield(event, 'status')
                            status = string(event.status);
                        end
                        app.updateSectionStatus(string(event.section), status);
                    end
                    if isfield(event,'role')&&strlength(string(event.role))>0
                        prefix(end+1)="role="+string(event.role); %#ok<AGROW>
                    end
                    if isfield(event,'phase')&&strlength(string(event.phase))>0
                        prefix(end+1)="phase="+string(event.phase); %#ok<AGROW>
                    end
                    if ~isempty(prefix)
                        message = "[" + strjoin(prefix, ", ") + "] " + message;
                    end
                else
                    message = string(event);
                end
                app.appendLog(message);
                drawnow limitrate;
            catch callback_error
                app.appendLog("[WARN] event callback: " + string(callback_error.message));
            end
        end

        function updateSectionStatus(app, section, status)
            app.StatusLabel.Text = char(section + " · " + status);
            data = app.SectionTable.Data;
            if ~istable(data) || ~ismember('Section', data.Properties.VariableNames)
                return;
            end
            idx = find(string(data.Section) == section, 1);
            if isempty(idx), return; end
            data.Status(idx) = status;
            app.SectionTable.Data = data;
        end

        function resetSectionStatuses(app)
            data = app.SectionTable.Data;
            if istable(data) && ismember('Status', data.Properties.VariableNames)
                data.Status(:) = "pending";
                if ismember('Action', data.Properties.VariableNames)
                    data.Status(string(data.Action) == "skip") = "skipped";
                end
                app.SectionTable.Data = data;
            end
        end

        function captureReceipt(app, receipt)
            if ~isstruct(receipt), return; end
            candidates = {'output_path','result_path','result_dir', ...
                'manifest_file','manifest_path'};
            for idx = 1:numel(candidates)
                name = candidates{idx};
                if isfield(receipt, name) && strlength(string(receipt.(name))) > 0
                    candidate = string(receipt.(name));
                    if ismember(name, {'manifest_file','manifest_path'})
                        candidate = string(fileparts(candidate));
                    end
                    if isfolder(candidate)
                        app.CurrentOutputPath = candidate;
                        break;
                    end
                end
            end
            if isfield(receipt, 'completed_sections')
                completed = string(receipt.completed_sections(:));
                for idx = 1:numel(completed)
                    app.updateSectionStatus(completed(idx), "completed");
                end
            end
            if isfield(receipt, 'failed_section') ...
                    && strlength(string(receipt.failed_section)) > 0
                app.updateSectionStatus(string(receipt.failed_section), "failed");
            end
            app.appendLog("Receipt: " + string(evalc('disp(receipt)')));
        end

        function setBusy(app, value)
            app.IsBusy = logical(value);
            if value
                app.StatusLamp.Color = [0.85 0.10 0.10];
            else
                app.StatusLamp.Color = [0.10 0.65 0.25];
            end
            app.RunButton.Enable = app.onOff(~value);
            app.ValidateButton.Enable = app.onOff(~value);
            app.LoadRequestButton.Enable = app.onOff(~value);
            app.ModeDropDown.Enable = app.onOff(~value);
            app.ScopeDropDown.Enable = app.onOff(~value);
            app.PresetDropDown.Enable = app.onOff(~value);
            if string(app.ScopeDropDown.Value) == "record"
                app.RecordAverageCheckBox.Enable = app.onOff(~value);
                app.RecordAverageOnlyCheckBox.Enable = app.onOff(~value && ...
                    app.RecordAverageCheckBox.Value);
            end
            app.SectionTable.Enable = app.onOff(~value);
            app.InputPathEdit.Enable = app.onOff(~value);
            app.BrowseInputButton.Enable = app.onOff(~value);
            app.ReuseSourceEdit.Enable = app.onOff(~value);
            app.BrowseReuseSourceButton.Enable = app.onOff(~value);
            app.DetectReuseSourceButton.Enable = app.onOff(~value);
            app.SaveRequestButton.Enable = app.onOff(~value);
            app.BasicSwitchTable.Enable = app.onOff(~value);
            app.BasicValueTable.Enable = app.onOff(~value);
            app.AdvancedSwitchTable.Enable = app.onOff(~value);
            app.AdvancedValueTable.Enable = app.onOff(~value);
            app.ExpertSwitchTable.Enable = app.onOff(~value);
            app.ExpertValueTable.Enable = app.onOff(~value);
            if value
                app.StatusLabel.Text = 'running';
            end
            drawnow;
        end

        function autoDetectScopeFromInput(app)
            path_value = string(app.InputPathEdit.Value);
            if ~isfolder(path_value), return; end
            [~,name] = fileparts(char(path_value));
            detected = "";
            if startsWith(string(name),"Cycle",'IgnoreCase',true)
                detected = "single";
            else
                try
                    entries = aaa.helpers.common.discover_record_cycles(path_value,[],true);
                    if ~isempty(entries), detected = "record"; end
                catch
                end
            end
            if strlength(detected)==0 || detected==string(app.ScopeDropDown.Value), return; end
            old_scope=string(app.Request.scope);
            previous=app.Request;
            mode = string(app.ModeDropDown.Value);
            next=aaa.schema.default_request(mode,detected);
            next.input_path=path_value;
            next.workflow=previous.workflow;
            next.params=previous.params;
            app.Request=next;
            app.ScopeDropDown.Value = detected;
            app.syncRequestToUI();
            app.appendLog("Input structure detected; scope changed: "+ ...
                old_scope+" -> "+detected+ ...
                ". Existing workflow and mode parameters were preserved.");
        end

        function quickInspectAndReport(app)
            app.appendLog("[QUICK CHECK] Inspecting input structure without result discovery.");
            drawnow limitrate;
            try
                plan=app.resolvedInspectionPlan();
                options=struct('scan_reuse',false);
                inspection=aaa.helpers.common.inspect_input_structure( ...
                    app.Request,plan,options);
                app.appendLog("[QUICK CHECK] path="+inspection.input_path+ ...
                    ", inferred_scope="+inspection.inferred_scope+ ...
                    ", cycles="+numel(inspection.cycles)+".");
                for idx=1:numel(inspection.cycles)
                    cycle=inspection.cycles(idx);
                    if cycle.empty,state="empty";elseif cycle.standard_input
                        state="movie input ready";else,state="movie input missing";end
                    app.appendLog("  "+cycle.label+": "+state+".");
                end
                for value=string(inspection.errors(:))'
                    app.appendLog("[ERROR] "+value);
                end
                app.appendLog("[QUICK CHECK] Reuse compatibility remains pending Detect/Validate/Run.");
            catch err
                app.appendLog("[ERROR] Quick input check failed: "+string(err.message));
            end
            drawnow limitrate;
        end

        function preflight=performAuthoritativePreflight(app,trigger)
            app.appendLog("[PREFLIGHT] Starting authoritative preflight: "+string(trigger)+".");
            drawnow limitrate;
            original_request=app.Request;
            preflight=struct('schema_version',"1.0.0", ...
                'created_at',datetime('now'),'trigger',string(trigger), ...
                'request_key',app.requestKey(original_request), ...
                'resolved_request_key',"",'valid',false, ...
                'inspection',struct(),'validation',struct(), ...
                'reuse_resolution',struct(),'resolved_request',original_request);
            try
                plan=app.resolvedInspectionPlan();
                progress=@(message)app.logPreflightProgress(message);
                options=struct('scan_reuse',true,'progress',progress);
                inspection=aaa.helpers.common.inspect_input_structure( ...
                    original_request,plan,options);
                app.reportInputInspection(inspection,plan);
                resolved_request=original_request;
                reuse_resolution=struct();
                if string(original_request.scope)=="single" && ...
                        ~isempty(inspection.cycles)
                    selected=string(inspection.cycles(1).selected_reuse_source);
                    if strlength(selected)>0
                        resolved_request.workflow.source_results_path=selected;
                        [resolved_request,reuse_resolution]= ...
                            aaa.helpers.common.resolve_request_reuse_source(resolved_request);
                    end
                end
                [resolved_request,validation]=aaa.schema.validate_request(resolved_request);
                app.reportValidation(validation,"权威 preflight");
                preflight.inspection=inspection;
                preflight.validation=validation;
                preflight.reuse_resolution=reuse_resolution;
                preflight.resolved_request=resolved_request;
                preflight.resolved_request_key=app.requestKey(resolved_request);
                preflight.valid=logical(inspection.valid)&&logical(validation.valid);
            catch err
                preflight.validation=struct('valid',false,'errors',string(err.message), ...
                    'warnings',strings(0,1));
                app.appendLog("[ERROR] Authoritative preflight failed: "+string(err.message));
            end
            app.AuthoritativePreflight=preflight;
            if preflight.valid,state="passed";else,state="failed";end
            app.appendLog("[PREFLIGHT] Authoritative record stored: "+state+ ...
                "; time="+string(preflight.created_at)+".");
            drawnow limitrate;
        end

        function [preflight,reused]=preflightForRun(app)
            reused=false;
            current_key=app.requestKey(app.Request);
            candidate=app.AuthoritativePreflight;
            if isstruct(candidate)&&isfield(candidate,'request_key')&& ...
                    string(candidate.request_key)==current_key
                if ~isfield(candidate,'valid')||~candidate.valid|| ...
                        app.preflightSourcesExist(candidate)
                    preflight=candidate;
                    reused=true;
                    return;
                end
                app.appendLog("[PREFLIGHT] A selected source path disappeared; validation record invalidated.");
                app.AuthoritativePreflight=struct();
            end
            preflight=app.performAuthoritativePreflight("Run fallback");
        end

        function tf=preflightSourcesExist(~,preflight)
            tf=true;
            if ~isfield(preflight,'inspection')||~isfield(preflight.inspection,'cycles'),return;end
            cycles=preflight.inspection.cycles;
            for idx=1:numel(cycles)
                source=string(cycles(idx).selected_reuse_source);
                if strlength(source)>0&&~(isfolder(source)||isfile(source))
                    tf=false;return;
                end
            end
        end

        function invalidatePreflight(app,reason)
            had_record=isstruct(app.AuthoritativePreflight)&& ...
                isfield(app.AuthoritativePreflight,'request_key');
            app.AuthoritativePreflight=struct();
            if had_record
                app.appendLog("[PREFLIGHT] Previous validation invalidated by "+string(reason)+".");
            else
                app.appendLog("[PREFLIGHT] Validation pending after "+string(reason)+".");
            end
        end

        function plan=resolvedInspectionPlan(app)
            [~,~,~,plan]=aaa.schema.resolve_sections(app.Request.mode, ...
                app.Request.workflow.preset,app.Request.workflow.sections);
            if app.Request.scope=="record"&&isfield(app.Request.record,'record_average_only') ...
                    &&app.Request.record.record_average_only
                [~,~,~,plan]=aaa.schema.resolve_sections( ...
                    app.Request.mode,"analysis_only",struct());
            end
        end

        function reportInputInspection(app,inspection,plan)
            app.appendLog("Input inspection: path="+inspection.input_path+ ...
                ", inferred_scope="+inspection.inferred_scope+ ...
                ", cycles="+numel(inspection.cycles)+".");
            plan_reuses=any(structfun(@(v)string(v)=="reuse",plan));
            for idx=1:numel(inspection.cycles)
                cycle=inspection.cycles(idx);
                if cycle.empty,state="empty";elseif cycle.standard_input
                    state="movie input ready";else,state="movie input missing";end
                app.appendLog("  "+cycle.label+": "+state+ ...
                    "; result folders="+cycle.output_candidate_count+ ...
                    ", compatible="+cycle.valid_output_count+".");
                if ~isempty(cycle.input_evidence)
                    app.appendLog("    input evidence: "+strjoin(cycle.input_evidence,", "));
                end
                if strlength(cycle.auto_reuse_candidate)>0
                    app.appendLog("    automatic candidate: "+cycle.auto_reuse_candidate);
                    app.logFrameRates("    detected saved frame rate",cycle.auto_frame_rates);
                end
                if strlength(inspection.manual_reuse_root)>0&& ...
                        strlength(cycle.manual_reuse_candidate)>0
                    app.appendLog("    manual candidate selected: "+ ...
                        cycle.manual_reuse_candidate+" [manual priority]");
                    app.logFrameRates("    manual saved frame rate",cycle.manual_frame_rates);
                elseif plan_reuses&&strlength(cycle.selected_reuse_source)>0
                    app.appendLog("    automatic Reuse source selected: "+ ...
                        cycle.selected_reuse_source);
                end
            end
            for value=string(inspection.warnings(:))',app.appendLog("[WARN] "+value);end
            for value=string(inspection.errors(:))',app.appendLog("[ERROR] "+value);end
            if inspection.valid
                app.appendLog("[OK] Input/reuse preflight passed.");
            else
                app.appendLog("[ERROR] Input/reuse preflight failed; Run is blocked.");
            end
        end

        function logPreflightProgress(app,message)
            app.appendLog("[PREFLIGHT] "+string(message));
            drawnow limitrate;
        end

        function text_value=preflightErrorText(app,preflight)
            lines=strings(0,1);
            if isfield(preflight,'inspection')&&isfield(preflight.inspection,'errors')
                lines=[lines;string(preflight.inspection.errors(:))]; %#ok<AGROW>
            end
            if isfield(preflight,'validation')&&isfield(preflight.validation,'errors')
                lines=[lines;string(preflight.validation.errors(:))]; %#ok<AGROW>
            end
            lines=unique(lines(strlength(lines)>0),'stable');
            if isempty(lines),lines="Authoritative preflight failed.";end
            text_value=char(strjoin(lines,newline));
        end

        function key=requestKey(~,request)
            try,key=string(jsonencode(request));catch,key="";end
        end

        function logButtonClick(app,action)
            app.appendLog("[UI] "+string(action)+" clicked.");
            drawnow limitrate;
        end

        function appendLog(app, message)
            message = string(message);
            message = splitlines(message);
            message = message(strlength(message) > 0);
            if isempty(message), return; end
            stamp = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
            app.LogLines = [app.LogLines; stamp + "  " + message]; %#ok<AGROW>
            if numel(app.LogLines) > 1200
                app.LogLines = app.LogLines(end-1199:end);
            end
            app.LogTextArea.Value = cellstr(app.LogLines);
            if ~isempty(app.SessionLogger) && isvalid(app.SessionLogger)
                level = "INFO";
                if any(startsWith(message, "[ERROR]"))
                    level = "ERROR";
                elseif any(startsWith(message, "[WARN]"))
                    level = "WARNING";
                end
                app.SessionLogger.write(level, strjoin(message, newline));
            end
            if isvalid(app.LogTextArea)
                try
                    scroll(app.LogTextArea, 'bottom');
                catch
                    % Scrolling support differs between MATLAB releases.
                end
            end
        end

        function logFrameRates(app,prefix,values)
            if ~isstruct(values),return;end
            roles=string(fieldnames(values));
            for rate_idx=1:numel(roles)
                value=values.(char(roles(rate_idx)));
                if isnumeric(value)&&isscalar(value)&&isfinite(value)
                    app.appendLog(prefix+": "+roles(rate_idx)+"="+value+" Hz.");
                end
            end
        end

        function UIFigureWindowButtonDown(app, ~)
            point = app.UIFigure.CurrentPoint;
            divider_position = app.LogDivider.Position;
            grid_position = app.MainGrid.Position;
            divider_bounds = [ ...
                grid_position(1) + divider_position(1), ...
                grid_position(2) + divider_position(2), ...
                divider_position(3), divider_position(4)];
            inside_divider = point(1) >= divider_bounds(1) ...
                && point(1) <= divider_bounds(1) + divider_bounds(3) ...
                && point(2) >= divider_bounds(2) ...
                && point(2) <= divider_bounds(2) + divider_bounds(4);
            if ~inside_divider, return; end
            app.LogResizeActive = true;
            app.LogResizeStartY = point(2);
            app.LogResizeStartHeight = app.LogPanel.Position(4);
            app.LogResizeAvailableHeight = app.MainTabs.Position(4) + ...
                app.LogPanel.Position(4);
            app.UIFigure.Pointer = 'top';
            app.UIFigure.WindowButtonMotionFcn = createCallbackFcn(app, ...
                @LogDividerMotion, true);
            app.UIFigure.WindowButtonUpFcn = createCallbackFcn(app, ...
                @LogDividerButtonUp, true);
        end

        function LogDividerMotion(app, ~)
            if ~app.LogResizeActive || ~isvalid(app.UIFigure), return; end
            point = app.UIFigure.CurrentPoint;
            requested = app.LogResizeStartHeight + ...
                (point(2) - app.LogResizeStartY);
            minimum_log_height = 160;
            minimum_tabs_height = 220;
            maximum_log_height = max(minimum_log_height, ...
                app.LogResizeAvailableHeight - minimum_tabs_height);
            log_height = min(max(requested,minimum_log_height), ...
                maximum_log_height);
            app.MainGrid.RowHeight = {205,'1x',22,log_height};
            drawnow limitrate nocallbacks;
        end

        function LogDividerButtonUp(app, ~)
            app.LogResizeActive = false;
            if isvalid(app.UIFigure)
                app.UIFigure.WindowButtonMotionFcn = [];
                app.UIFigure.WindowButtonUpFcn = [];
                app.UIFigure.Pointer = 'arrow';
            end
        end

        function value = parseTextValue(app, raw, entry)
            raw = strtrim(string(raw));
            value_type = string(entry.value_type);
            switch value_type
                case {"string", "path"}
                    value = raw;
                case "string_vector"
                    if strlength(raw) == 0 || raw == "[]"
                        value = strings(0, 1);
                    else
                        cleaned = regexprep(raw, '[,;]', ' ');
                        value = split(strtrim(cleaned));
                        value = value(strlength(value) > 0);
                    end
                case "enum"
                    if ~ismember(raw, entry.choices)
                        error('AAA:App:InvalidEnum', '允许值: %s', ...
                            strjoin(entry.choices, ', '));
                    end
                    value = raw;
                case {"double", "integer"}
                    if strlength(raw) == 0 || raw == "[]"
                        value = [];
                    else
                        value = str2double(raw);
                        if isnan(value) && lower(raw) ~= "nan"
                            error('AAA:App:InvalidNumber', '请输入数值、NaN、Inf 或 []。');
                        end
                        if value_type == "integer" && isfinite(value) && value ~= fix(value)
                            error('AAA:App:InvalidInteger', '该参数必须为整数。');
                        end
                    end
                case "numeric_vector"
                    value = app.parseNumericVector(raw);
                case "struct"
                    if strlength(raw) == 0 || ismember(lower(raw), ["struct()", "{}"])
                        value = struct();
                    else
                        value = jsondecode(char(raw));
                        if ~isstruct(value)
                            error('AAA:App:InvalidStruct', ...
                                '结构体参数请输入 JSON object，例如 {""field"":1}。');
                        end
                    end
                otherwise
                    value = raw;
            end
            if isnumeric(value) && ~isempty(value)
                if ~isempty(entry.minimum) && any(value < entry.minimum, 'all')
                    error('AAA:App:BelowMinimum', '最小值为 %s。', mat2str(entry.minimum));
                end
                if ~isempty(entry.maximum) && any(value > entry.maximum, 'all')
                    error('AAA:App:AboveMaximum', '最大值为 %s。', mat2str(entry.maximum));
                end
            end
        end

        function values = parseNumericVector(~, raw)
            if strlength(raw) == 0 || raw == "[]"
                values = [];
                return;
            end
            cleaned = regexprep(raw, '[\[\],;]', ' ');
            tokens = split(strtrim(cleaned));
            tokens = tokens(strlength(tokens) > 0);
            values = zeros(1, numel(tokens));
            for idx = 1:numel(tokens)
                values(idx) = str2double(tokens(idx));
                if isnan(values(idx)) && lower(tokens(idx)) ~= "nan"
                    error('AAA:App:InvalidVector', ...
                        '向量必须类似 [1 2 3]、[0 0] 或 []。');
                end
            end
        end

        function text_value = valueToText(~, value)
            if isempty(value)
                text_value = "[]";
            elseif isstring(value)
                if isscalar(value), text_value = value; else, text_value = strjoin(value, " "); end
            elseif ischar(value)
                text_value = string(value);
            elseif isnumeric(value) || islogical(value)
                text_value = string(mat2str(value));
            elseif isstruct(value)
                if isempty(fieldnames(value))
                    text_value = "{}";
                else
                    text_value = string(jsonencode(value));
                end
            else
                text_value = string(evalc('disp(value)'));
            end
        end

        function UIFigureCloseRequest(app, ~)
            if app.IsBusy
                choice = uiconfirm(app.UIFigure, ...
                    ['分析正在运行。关闭界面不会安全地中止正在执行的科学函数，' ...
                     '因此当前不允许关闭。'], ...
                    'AAA 正在运行', 'Options', {'继续等待'}, ...
                    'DefaultOption', 1, 'CancelOption', 1);
                if ~isempty(choice), return; end
            end
            delete(app);
        end

        function createParameterTab(app, tab, level)
            grid = uigridlayout(tab, [4 1]);
            grid.RowHeight = {24, '1x', 24, '2x'};
            grid.Padding = [8 8 8 8];
            switch_label = uilabel(grid, 'Text', ...
                '公开判断与默认开关（Enabled 列为真实复选框）');
            switch_table = uitable(grid);
            switch_table.ColumnName = {'参数','Enabled','分组','Section','Key','说明/条件/来源'};
            switch_table.ColumnEditable = [false true false false false false];
            switch_table.ColumnWidth = {190,70,120,110,240,'auto'};
            switch_table.UserData = level + "-switch";
            switch_table.CellEditCallback = ...
                @(src,event) SwitchTableCellEdit(app,src,event);

            value_label = uilabel(grid, 'Text', ...
                '公开数值、路径与选项（Value 列可编辑，Allowed 显示枚举值）');
            value_table = uitable(grid);
            value_table.ColumnName = {'参数','Value','单位','Allowed','分组','Section','Key','说明/条件/来源'};
            value_table.ColumnEditable = [false true false false false false false false];
            value_table.ColumnWidth = {185,145,70,180,110,100,235,'auto'};
            value_table.UserData = level + "-value";
            value_table.CellEditCallback = ...
                @(src,event) ValueTableCellEdit(app,src,event);

            switch level
                case "Basic"
                    app.BasicGrid = grid;
                    app.BasicSwitchLabel = switch_label;
                    app.BasicSwitchTable = switch_table;
                    app.BasicValueLabel = value_label;
                    app.BasicValueTable = value_table;
                case "Advanced"
                    app.AdvancedGrid = grid;
                    app.AdvancedSwitchLabel = switch_label;
                    app.AdvancedSwitchTable = switch_table;
                    app.AdvancedValueLabel = value_label;
                    app.AdvancedValueTable = value_table;
                case "Expert"
                    app.ExpertGrid = grid;
                    app.ExpertSwitchLabel = switch_label;
                    app.ExpertSwitchTable = switch_table;
                    app.ExpertValueLabel = value_label;
                    app.ExpertValueTable = value_table;
            end
        end

        function createComponents(app)
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [40 25 1500 950];
            app.UIFigure.Name = 'AP Analysis APP (AAA)';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, ...
                @UIFigureCloseRequest, true);
            app.UIFigure.WindowButtonDownFcn = createCallbackFcn(app, ...
                @UIFigureWindowButtonDown, true);

            app.MainGrid = uigridlayout(app.UIFigure, [4 1]);
            app.MainGrid.RowHeight = {225, '1x', 22, 320};
            app.MainGrid.Padding = [10 10 10 10];
            app.MainGrid.RowSpacing = 6;

            app.HeaderPanel = uipanel(app.MainGrid, 'Title', ...
                'AAA · request 控制层（算法结果保存在目标目录）');
            app.HeaderPanel.Layout.Row = 1;
            app.HeaderGrid = uigridlayout(app.HeaderPanel, [4 10]);
            app.HeaderGrid.RowHeight = {30, 30, 50, 30};
            app.HeaderGrid.ColumnWidth = {90,110,70,110,70,'1x',80,105,105,105};
            app.HeaderGrid.Padding = [10 10 10 10];
            app.HeaderGrid.RowSpacing = 5;

            app.TitleLabel = uilabel(app.HeaderGrid, 'Text', ...
                'AP Analysis APP · AAA');
            app.TitleLabel.FontSize = 18;
            app.TitleLabel.FontWeight = 'bold';
            app.TitleLabel.Layout.Row = 1;
            app.TitleLabel.Layout.Column = [1 6];
            app.StatusLamp = uilamp(app.HeaderGrid);
            app.StatusLamp.Color = [0.10 0.65 0.25];
            app.StatusLamp.Tooltip = 'Green: idle; Red: analysis running';
            app.StatusLamp.Layout.Row = 1;
            app.StatusLamp.Layout.Column = 7;
            app.StatusLabel = uilabel(app.HeaderGrid, 'Text', 'starting');
            app.StatusLabel.HorizontalAlignment = 'right';
            app.StatusLabel.FontWeight = 'bold';
            app.StatusLabel.Layout.Row = 1;
            app.StatusLabel.Layout.Column = [8 10];

            app.ModeLabel = uilabel(app.HeaderGrid, 'Text', '模式');
            app.ModeLabel.Layout.Row = 2; app.ModeLabel.Layout.Column = 1;
            app.ModeDropDown = uidropdown(app.HeaderGrid);
            app.ModeDropDown.Items = {'AP','Dual'};
            app.ModeDropDown.ItemsData = ["ap","dual"];
            app.ModeDropDown.Layout.Row = 2; app.ModeDropDown.Layout.Column = 2;
            app.ModeDropDown.ValueChangedFcn = createCallbackFcn(app, ...
                @ModeDropDownValueChanged, true);

            app.ScopeLabel = uilabel(app.HeaderGrid, 'Text', '范围');
            app.ScopeLabel.Layout.Row = 2; app.ScopeLabel.Layout.Column = 3;
            app.ScopeDropDown = uidropdown(app.HeaderGrid);
            app.ScopeDropDown.Items = {'Single cycle','Record'};
            app.ScopeDropDown.ItemsData = ["single","record"];
            app.ScopeDropDown.Layout.Row = 2; app.ScopeDropDown.Layout.Column = 4;
            app.ScopeDropDown.ValueChangedFcn = createCallbackFcn(app, ...
                @ScopeDropDownValueChanged, true);

            app.InputLabel = uilabel(app.HeaderGrid, 'Text', '输入目录');
            app.InputLabel.Layout.Row = 2; app.InputLabel.Layout.Column = 5;
            app.InputPathEdit = uieditfield(app.HeaderGrid, 'text');
            app.InputPathEdit.Placeholder = 'Cycle 或 Record 路径';
            app.InputPathEdit.Layout.Row = 2; app.InputPathEdit.Layout.Column = 6;
            app.InputPathEdit.ValueChangedFcn = createCallbackFcn(app, ...
                @InputPathEditValueChanged, true);
            app.BrowseInputButton = uibutton(app.HeaderGrid, 'push', 'Text', '浏览…');
            app.BrowseInputButton.Layout.Row = 2; app.BrowseInputButton.Layout.Column = 7;
            app.BrowseInputButton.ButtonPushedFcn = createCallbackFcn(app, ...
                @BrowseInputButtonPushed, true);

            app.LoadRequestButton = uibutton(app.HeaderGrid, 'push', 'Text', '载入 request');
            app.LoadRequestButton.Layout.Row = 2; app.LoadRequestButton.Layout.Column = 8;
            app.LoadRequestButton.ButtonPushedFcn = createCallbackFcn(app, ...
                @LoadRequestButtonPushed, true);
            app.SaveRequestButton = uibutton(app.HeaderGrid, 'push', 'Text', '保存 request');
            app.SaveRequestButton.Layout.Row = 2; app.SaveRequestButton.Layout.Column = 9;
            app.SaveRequestButton.ButtonPushedFcn = createCallbackFcn(app, ...
                @SaveRequestButtonPushed, true);
            app.OpenOutputButton = uibutton(app.HeaderGrid, 'push', 'Text', '打开输出');
            app.OpenOutputButton.Layout.Row = 2; app.OpenOutputButton.Layout.Column = 10;
            app.OpenOutputButton.ButtonPushedFcn = createCallbackFcn(app, ...
                @OpenOutputButtonPushed, true);

            app.PresetLabel = uilabel(app.HeaderGrid, 'Text', 'Preset');
            app.PresetLabel.Layout.Row = 3; app.PresetLabel.Layout.Column = 1;
            app.PresetDropDown = uidropdown(app.HeaderGrid);
            app.PresetDropDown.Items = {'Full','Reuse motion','Reuse trace', ...
                'Skip motion','Analysis only','Motion only','Custom'};
            app.PresetDropDown.ItemsData = ["full","reuse_motion","reuse_trace", ...
                "skip_motion","analysis_only","motion_only","custom"];
            app.PresetDropDown.Layout.Row = 3; app.PresetDropDown.Layout.Column = [2 4];
            app.PresetDropDown.ValueChangedFcn = createCallbackFcn(app, ...
                @PresetDropDownValueChanged, true);

            app.RecordAnalysisPanel = uipanel(app.HeaderGrid, ...
                'Title','Record analysis');
            app.RecordAnalysisPanel.Layout.Row = 3;
            app.RecordAnalysisPanel.Layout.Column = [5 7];
            app.RecordAnalysisGrid = uigridlayout(app.RecordAnalysisPanel,[1 2]);
            app.RecordAnalysisGrid.ColumnWidth = {'1x','1x'};
            app.RecordAnalysisGrid.Padding = [6 2 6 2];
            app.RecordAverageCheckBox = uicheckbox(app.RecordAnalysisGrid, ...
                'Text','Run Record-level analysis');
            app.RecordAverageCheckBox.ValueChangedFcn = createCallbackFcn(app, ...
                @RecordAverageCheckBoxValueChanged,true);
            app.RecordAverageOnlyCheckBox = uicheckbox(app.RecordAnalysisGrid, ...
                'Text','Use existing Cycle results only');
            app.RecordAverageOnlyCheckBox.ValueChangedFcn = createCallbackFcn(app, ...
                @RecordAverageOnlyCheckBoxValueChanged,true);

            app.ValidateButton = uibutton(app.HeaderGrid, 'push', 'Text', '验证 request');
            app.ValidateButton.Layout.Row = 3; app.ValidateButton.Layout.Column = 8;
            app.ValidateButton.ButtonPushedFcn = createCallbackFcn(app, ...
                @ValidateButtonPushed, true);
            app.RunButton = uibutton(app.HeaderGrid, 'push', 'Text', '运行 AAA');
            app.RunButton.FontWeight = 'bold';
            app.RunButton.BackgroundColor = [0.20 0.55 0.30];
            app.RunButton.FontColor = [1 1 1];
            app.RunButton.Layout.Row = 3; app.RunButton.Layout.Column = [9 10];
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, ...
                @RunButtonPushed, true);

            app.ReuseSourceLabel = uilabel(app.HeaderGrid, 'Text', 'Reuse source');
            app.ReuseSourceLabel.Layout.Row = 4;
            app.ReuseSourceLabel.Layout.Column = 1;
            app.ReuseSourceEdit = uieditfield(app.HeaderGrid, 'text');
            app.ReuseSourceEdit.Placeholder = ...
                'Result folder, Cycle/Record search root, or blank for automatic detection';
            app.ReuseSourceEdit.Layout.Row = 4;
            app.ReuseSourceEdit.Layout.Column = [2 6];
            app.ReuseSourceEdit.ValueChangedFcn = createCallbackFcn(app, ...
                @ReuseSourceEditValueChanged, true);
            app.BrowseReuseSourceButton = uibutton(app.HeaderGrid, 'push', ...
                'Text', 'Browse reuse');
            app.BrowseReuseSourceButton.Layout.Row = 4;
            app.BrowseReuseSourceButton.Layout.Column = 7;
            app.BrowseReuseSourceButton.ButtonPushedFcn = createCallbackFcn(app, ...
                @BrowseReuseSourceButtonPushed, true);
            app.DetectReuseSourceButton = uibutton(app.HeaderGrid, 'push', ...
                'Text', 'Detect reuse');
            app.DetectReuseSourceButton.Layout.Row = 4;
            app.DetectReuseSourceButton.Layout.Column = 8;
            app.DetectReuseSourceButton.ButtonPushedFcn = createCallbackFcn(app, ...
                @DetectReuseSourceButtonPushed, true);

            app.MainTabs = uitabgroup(app.MainGrid);
            app.MainTabs.Layout.Row = 2;
            app.WorkflowTab = uitab(app.MainTabs, 'Title', 'Sections');
            app.WorkflowGrid = uigridlayout(app.WorkflowTab, [2 1]);
            app.WorkflowGrid.RowHeight = {32, '1x'};
            app.WorkflowGrid.Padding = [8 8 8 8];
            app.WorkflowHelpLabel = uilabel(app.WorkflowGrid, 'Text', ...
                ['Preset 只初始化计划；手工修改 Action 后 request 变为 custom。' ...
                 ' Status 由真实 section/cycle 事件更新。']);
            app.SectionTable = uitable(app.WorkflowGrid);
            app.SectionTable.ColumnName = {'Section','名称','Action','Status','说明'};
            app.SectionTable.ColumnEditable = [false false true false false];
            app.SectionTable.ColumnWidth = {130,170,100,110,'auto'};
            app.SectionTable.CellEditCallback = createCallbackFcn(app, ...
                @SectionTableCellEdit, true);

            app.BasicTab = uitab(app.MainTabs, 'Title', 'Basic');
            app.createParameterTab(app.BasicTab, "Basic");
            app.AdvancedTab = uitab(app.MainTabs, 'Title', 'Advanced');
            app.createParameterTab(app.AdvancedTab, "Advanced");
            app.ExpertTab = uitab(app.MainTabs, 'Title', 'Expert');
            app.createParameterTab(app.ExpertTab, "Expert");

            app.LogDivider = uilabel(app.MainGrid, ...
                'Text','↕ 拖动调整 Sections / Run Log 高度');
            app.LogDivider.Layout.Row = 3;
            app.LogDivider.HorizontalAlignment = 'center';
            app.LogDivider.FontColor = [0.25 0.30 0.36];
            app.LogDivider.BackgroundColor = [0.88 0.90 0.93];
            app.LogDivider.Tooltip = '按住并上下拖动';
            app.LogPanel = uipanel(app.MainGrid,'Title','Run Log（可拖动调整）');
            app.LogPanel.Layout.Row = 4;
            app.LogGrid = uigridlayout(app.LogPanel, [2 1]);
            app.LogGrid.RowHeight = {28, '1x'};
            app.LogGrid.Padding = [8 6 8 8];
            app.LogHelpLabel = uilabel(app.LogGrid, 'Text', ...
                ['这里只显示 App 会话和后端事件；目标目录中的 analysis.log ' ...
                 '与 analysis_manifest.mat 才是正式追溯记录。']);
            app.LogTextArea = uitextarea(app.LogGrid, 'Editable', 'off');
            app.LogTextArea.FontName = 'Consolas';

            app.UIFigure.Visible = 'on';
        end
    end

    methods (Static, Access = private)
        function [value, found] = getPath(input, path_text)
            parts = cellstr(split(string(path_text), '.'));
            value = input;
            found = true;
            for idx = 1:numel(parts)
                name = parts{idx};
                if ~isstruct(value) || ~isscalar(value) || ~isfield(value, name)
                    value = [];
                    found = false;
                    return;
                end
                value = value.(name);
            end
        end

        function output = setPath(input, path_text, value)
            parts = cellstr(split(string(path_text), '.'));
            output = AP_Analysis_APP.setNested(input, parts, value);
        end

        function output = setNested(input, parts, value)
            output = input;
            name = parts{1};
            if numel(parts) == 1
                output.(name) = value;
                return;
            end
            if ~isfield(output, name) || ~isstruct(output.(name)) || ...
                    ~isscalar(output.(name))
                output.(name) = struct();
            end
            output.(name) = AP_Analysis_APP.setNested( ...
                output.(name), parts(2:end), value);
        end

        function value = onOff(flag)
            if flag, value = 'on'; else, value = 'off'; end
        end
    end

    methods (Access = public)
        function app = AP_Analysis_APP
            createComponents(app);
            registerApp(app, app.UIFigure);
            runStartupFcn(app, @startupFcn);
            if nargout == 0
                clear app
            end
        end

        function delete(app)
            if ~isempty(app.SessionLogger) && isvalid(app.SessionLogger)
                app.SessionLogger.close();
            end
            if ~isempty(app.UIFigure) && isvalid(app.UIFigure)
                delete(app.UIFigure);
            end
        end
    end
end
