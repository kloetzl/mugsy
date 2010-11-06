' SEQAN Development Build System
' Visual Studio Macro Module
'
' 2003/2004/2008 by Andreas Doering
' doering@inf.fu-berlin.de
'________________________________________________________________________
'
' Installing:
' 1. Import this macro file into your Visual Studio (version 7.10 +)
' 2. Write Macros that use CallMacro(.) to execute commands 
'    (or adjust CallMacroGCC() and CallMacroSUN())
' 3. Customize the toolbar of Visual Studio for a quick macro access
' 4. Execute UpdateProjectView() to update the project view and configurations
'________________________________________________________________________
'
' Macros:
' - UpdateProjectView():                update the project view
' - EnterPassWord():                    display password dialog
' - CallMacro(command, need_passwd):    execute command
'
' Helper Functions:
' - CreatePlinkCommand(platform, server, path, tool): create plink command string
'________________________________________________________________________
'

Imports EnvDTE
Imports System
Imports System.Diagnostics
Imports System.Windows.Forms
Imports System.Collections
Imports Microsoft.VisualBasic.FileSystem


Public Module Seqan
    ' user configuration
    Dim g_HostGCC As String = "fichte"
    Dim g_HostSUN As String = "human"
    Dim g_SeqanRootGCC As String = "~/seqan/"
    Dim g_SeqanRootSUN As String = g_SeqanRootGCC
    Dim g_PythonBin As String = "C:\Program Files\Python\python.exe"

    Dim g_PassWord
    Dim g_PassWordForm As New Form
    Dim g_PassWordTextBox As New TextBox
    Dim g_PassWordButtonOK As New Button
    Dim g_IgnoreDirs As String = "/.svn/CVS/windows/gcc/sun/mingw/lara/"

    Sub CallDDDOC()
        Dim solutionpath As String = DTE.Solution.Properties.Item("Path").Value
        Dim pos As Integer = InStrRev(solutionpath, "\")
        If (pos) Then
            solutionpath = Mid(solutionpath, 1, pos - 1)
        End If
        Dim dddocpath As String = solutionpath & "\..\docs\dddoc"
        Dim librarypath As String = solutionpath & "\..\projects\library"

        ExecuteCommand("Tools.Shell /o /dir:""" & dddocpath & """ " & g_PythonBin & " " & dddocpath & "\main.py """ & librarypath & """")
        'ExecuteCommand("Tools.Shell /o /dir:""" & dddocpath & """ cmd /c " & dddocpath & "\main.py """ & librarypath & """")
    End Sub
    ' GCC Macro
    ' todo: adjust command dtring
    Sub CallMacroGCC()
        CallMacro(CreatePlinkCommand("gcc", g_HostGCC, g_SeqanRootGCC, "make", "build"), True)
    End Sub

    ' GCC Macro
    ' todo: adjust command dtring
    Sub CallMacroGCCRun()
        CallMacro(CreatePlinkCommand("gcc", g_HostGCC, g_SeqanRootGCC, "make", "run"), True)
    End Sub

    ' Make Clean Command
    Sub CallMacroClean()
        CallMacro(CreatePlinkCommand("gcc", g_HostGCC, g_SeqanRootGCC, "make", "clean"), True)
    End Sub

    ' SUN Macro
    ' todo: adjust command dtring
    Sub CallMacroSUN()
        CallMacro(CreatePlinkCommand("sun", g_HostSUN, g_SeqanRootSUN, "/import/gnu/bin/make", "run"), True)
    End Sub

    'Create PLINK Command String
    Private Function CreatePlinkCommand(ByVal platform As String, ByVal server As String, ByVal path As String, ByVal tool As String, ByVal command As String)
        Return "plink.exe -batch -ssh -l %username% -pw %password% " & server & " " & tool & " -s -C """ & path & """ " & command & " Platform=" & platform & " Project=%project% Mode=%mode%"
    End Function

    'Execute Command
    Private Sub CallMacro(ByVal command As String, ByVal need_passwd As Boolean)
        Dim project_mode As String = DTE.Solution.SolutionBuild.ActiveConfiguration.Name
        Dim solutionpath As String = DTE.Solution.Properties.Item("Path").Value
        Dim pos As Integer = InStrRev(solutionpath, "\")
        If (pos) Then
            solutionpath = Mid(solutionpath, 1, pos - 1)
        End If

        If (need_passwd And (g_PassWord = "")) Then
            EnterPassWord()
        End If

        If (Not need_passwd Or (g_PassWord <> "")) Then
            ExecuteCommand("Tools.Shell /o /dir:""" & solutionpath & """ CScript //Nologo Macro.wsf  """ & command & """ " & project_mode & " """ & g_PassWord & """")
        End If
    End Sub


    ' creates new configuration
    Private Sub CreateConfiguration(ByVal path As String, ByVal name As String, ByVal mode As String)
        On Error Resume Next

        'create solution configuration
        Dim solconfigs As SolutionConfigurations = DTE.Solution.SolutionBuild.SolutionConfigurations
        Dim solconfig As SolutionConfiguration = solconfigs.Add(name & " " & mode, mode, True)
        solconfig.Activate()

        'set project configuration properties
        Dim proj As Project = DTE.ActiveSolutionProjects(0)
        Dim config As EnvDTE.Configuration = proj.ConfigurationManager.ActiveConfiguration
        config.Properties.Item("OutputDirectory").Value = name
        config.Properties.Item("IntermediateDirectory").Value = path & "windows\"
        config.Properties.Item("WorkingDirectory").Value = ".."
    End Sub

    Private Function IsExpandFolder(ByVal level, ByVal file)
        Return ((level <= 2) Or (InStr(g_IgnoreDirs, "/" & file & "/") = 0))
    End Function

    ' Recursive Update Project View
    Private Sub UpdateProjectView_1(ByRef item As ProjectItem, ByVal path As String, ByVal level As Integer, ByVal is_apps_subfolder As Boolean)
        'create configuration (if item is a project)
        If (((level = 3) And (item.Name <> "demos") And (item.Name <> "apps") And (item.Name <> "seqan")) Or ((level = 4) And is_apps_subfolder)) Then
            'new_configs.Add(item.Name, path)
            CreateConfiguration(path, item.Name, "Debug")
            CreateConfiguration(path, item.Name, "Release")
        End If

        'search for files
        Dim file
        Dim file_list As New ArrayList
        Dim folder_list As New ArrayList
        file = Dir(path, FileAttribute.Directory)
        Do While (file <> "")
            If (file <> ".svn") And (file <> "CVS") Then
                If (GetAttr(path & file) And FileAttribute.Directory) Then
                    If (IsExpandFolder(level, file)) Then
                        folder_list.Add(file)
                    End If
                Else
                    If file.IndexOf(".h") > 0 Or file.IndexOf(".cpp") > 0 Then
                        file_list.Add(file)
                    End If
                End If
            End If
            file = Dir()
        Loop

        'insert files into treeview
        file_list.Sort()
        For i As Integer = 0 To file_list.Count - 1
            Try
                item.ProjectItems.AddFromFile(path & file_list(i))
            Catch ex As Exception
                'MsgBox("file: " & path & file_list(i))
            End Try
        Next

        'insert folder into treeview
        folder_list.Sort()
        For i As Integer = 0 To folder_list.Count - 1
            Try
                Dim new_folder As ProjectItem = item.ProjectItems.AddFolder(folder_list(i), "{6BB5F8F0-4483-11D3-8BCF-00C04F8EC28C}")
                UpdateProjectView_1(new_folder, path & folder_list(i) & "\", level + 1, item.Name = "apps")
            Catch ex As Exception
                'MsgBox("folder: " & folder_list(i))
            End Try
        Next

    End Sub

    ' Updates Project View and List of Configurations
    Sub UpdateProjectView()
        On Error Resume Next

        Dim proj As Project = DTE.ActiveSolutionProjects(0)
        Dim configs As SolutionConfigurations = DTE.Solution.SolutionBuild.SolutionConfigurations

        'stores active configuration
        Dim activconfigname As String = DTE.Solution.SolutionBuild.ActiveConfiguration.Name

        'empty configuration list
        Dim config As SolutionConfiguration
        For Each config In configs
            If (config.Name <> "Debug") And (config.Name <> "Release") Then
                proj.ConfigurationManager.DeleteConfigurationRow(config.Name)
                config.Delete()
            End If
        Next

        'empty project view
        Dim item As ProjectItem
        For Each item In proj.ProjectItems
            If (item.Name = "projects") Then
                item.Remove()
            End If
        Next

        'create root
        item = proj.ProjectItems.AddFolder("projects", "{6BB5F8F0-4483-11D3-8BCF-00C04F8EC28C}")

        'get root path
        Dim solutionpath As String = DTE.Solution.Properties.Item("Path").Value
        Dim pos As Integer = InStrRev(solutionpath, "\")
        If (pos) Then
            solutionpath = Mid(solutionpath, 1, pos - 1)
        End If
        pos = InStrRev(solutionpath, "\")
        If (pos) Then
            solutionpath = Mid(solutionpath, 1, pos - 1)
        End If

        'start recursion to built up tree
        UpdateProjectView_1(item, solutionpath & "\projects\", 1, False)

        'restore active solution
        DTE.Solution.SolutionBuild.SolutionConfigurations.Item(activconfigname).Activate()
    End Sub

    ' Handler for OK Button of Password Dialog
    Private Sub g_PassWordButtonOK_Click(ByVal sender As System.Object, ByVal e As System.EventArgs)
        g_PassWord = g_PassWordTextBox.Text
        g_PassWordForm.Close()
    End Sub


    ' Display Password Dialog
    Sub EnterPassWord()

        g_PassWordForm.FormBorderStyle = FormBorderStyle.FixedToolWindow
        g_PassWordForm.Text = "RBuildTool - Enter Password"
        g_PassWordForm.Width = 225
        g_PassWordForm.Height = 95
        g_PassWordForm.TopMost = True

        g_PassWordButtonOK.Text = "OK"
        g_PassWordButtonOK.Top = 40
        g_PassWordButtonOK.Left = 135

        g_PassWordTextBox.Width = 200
        g_PassWordTextBox.Left = 10
        g_PassWordTextBox.Top = 10
        g_PassWordTextBox.PasswordChar = "*"

        g_PassWordForm.Controls.Add(g_PassWordButtonOK)
        g_PassWordForm.Controls.Add(g_PassWordTextBox)

        g_PassWordForm.AcceptButton = g_PassWordButtonOK
        AddHandler g_PassWordButtonOK.Click, New System.EventHandler(AddressOf g_PassWordButtonOK_Click)

        g_PassWordForm.ShowDialog()

    End Sub



End Module
