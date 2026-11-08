#!/usr/bin/env perl

# Provisional end-to-end launcher for AbRaCD.
use strict;
use warnings;
use FindBin qw($Bin);
use File::Basename qw(basename dirname);
use File::Spec;
use Cwd qw(getcwd);
use POSIX qw(strftime);
use IO::Handle;
use Fcntl qw(O_WRONLY O_CREAT O_EXCL);
use Errno qw(EEXIST);

my $PREPROCESS = File::Spec->catfile($Bin, 'AbRaCD_chothia_preprocces.pl');
my $ENGINEERING_PREP = File::Spec->catfile($Bin, 'ab_engineering_prep.pl');
my $REFINE     = File::Spec->catfile($Bin, 'refine_loops.pl');
my $REPO_ROOT  = File::Spec->rel2abs(File::Spec->catdir($Bin, '..', '..'));
my $MAPS_DIR   = File::Spec->catdir($REPO_ROOT, 'data', 'Maps');
my $CDR_MAP    = File::Spec->catfile($MAPS_DIR, 'CDRS.bin');
my $KORP_MAP   = File::Spec->catfile($MAPS_DIR, 'korp6Dv1.bin');
my $DUN_MAP    = File::Spec->catfile($MAPS_DIR, 'dunbrack.bin');

# Default RCD settings, based on scripts/bash/script_allshort_CDRs_RCD.sh.
my $NP       = 1;
my %CDR_ARGS = (
    H1 => ['--H1', 0.5, 0.2, 20000, 2000], H2 => ['--H2', 0.9, 0.5, 20000, 2000],
    H3 => ['--H3', 0.5, 0.5, 50000, 1000], L1 => ['--L1', 0.9, 0.5, 20000, 2000],
    L2 => ['--L2', 0.9, 0.5, 20000, 2000], L3 => ['--L3', 0.9, 0.5, 20000, 2000],
);
my %CDR_EXTRA_ARGS = (
    H3 => ['--kink_merge', 0.9, '--kink_mergeA', 0.9, '--tkink', 0.70, '--tkinkA', 0.70],
);

our ($LOG_FILE, $LOG_FH, $CURRENT_STAGE, $RERUN_COMMAND, @TEMP_RCD_INPUT_FILES);
$CURRENT_STAGE = 'argument validation';

sub quote_arg {
    my ($value) = @_;
    $value = '' unless defined $value;
    $value =~ s/'/'\\''/g;
    return "'$value'";
}

sub format_command {
    return join(' ', map { quote_arg($_) } @_);
}

sub log_event {
    my ($level, $message) = @_;
    return unless $LOG_FH;
    my $timestamp = strftime('%Y-%m-%d %H:%M:%S', localtime);
    $message =~ s/\n/\n                      /g;
    print {$LOG_FH} "[$timestamp] [$level] [$CURRENT_STAGE] $message\n";
}

sub cleanup_temp_rcd_inputs {
    for my $path (@TEMP_RCD_INPUT_FILES) {
        next unless defined($path) && -e $path;
        if (unlink $path) {
            log_event('CLEANUP', "Removed temporary RCD input: $path");
        } else {
            log_event('WARNING', "Could not remove temporary RCD input '$path': $!");
        }
    }
    @TEMP_RCD_INPUT_FILES = ();
}

sub usage {
    my ($fh) = @_;
    print {$fh} <<'USAGE';
abracd.pl — run the AbRaCD workflow for one or more CDRs

Usage:
  abracd.pl [STRUCTURE.pdb | LIST.txt] [--raw | --model] [--purge] --CDRs CDR [CDR ... | ALL] [OPTIONS]
  abracd.pl [--raw | --model] --CDRs CDR [CDR ... | ALL] [OPTIONS]

Input:
  No structure/list:   process all PDB files in the current directory.
  STRUCTURE.pdb:       process one structure.
  LIST.txt:            one PDB per line (relative paths use the list directory).

  --CDRs:              required. Provide one or more of H1, H2, H3, L1, L2, L3,
                       or provide ALL by itself. The first CDR listed is processed
                       first. Without --raw, preprocessing runs once before the
                       first CDR and later CDRs reuse its generated lists.
                       ALL expands to H1 H2 H3 L1 L2 L3.

  --purge:             use the same input file as the purge file. It must have
                       one row per structure: PDB  H_chain  L_chain (use - for
                       a missing chain). Only applies without --raw.

  --raw:               skip Chothia preprocessing; place immediately after
                       the input argument. With one CDR, the input argument is
                       required and must be its RCD input list. With multiple
                       CDRs or ALL, do not provide an input argument; the script
                       finds the required AbRaCD_input_<CDR>.txt files in the
                       current directory.

  --OUT NAME:          optional output name prefix (default: none). Output folders
                       and Rosetta suffixes include NAME when provided.

  --model:             prepare and run an antibody design/modeling workflow
                       without a known native structure. Requires one RCD input
                       list and exactly one CDR. Runs engineering preprocessing,
                       RCD, and Rosetta refinement in modeling mode, then removes
                       native-dependent RMSD columns. Do not combine with --raw
                       or --purge.

  --NP X:              use the MPI RCD executable with X processes (default: serial RCD).

  --NLOOPS X:          RCD loop sampling count (default: 50000 for H3; 20000 otherwise).

  --LOCOBEST X:        RCD loco_best count (default: 1000 for H3; 2000 otherwise).

  --repRCD X:          number of RCD runs (default: 1).

  --repROS Y:          Rosetta refinement repetitions per RCD run (default: 1; 0 skips refinement).
  
  Log file:            NAME.log when --OUT NAME is used; otherwise AbRaCD.log.
                       The log is written in the current working directory.

Examples:
  abracd.pl antibody.pdb --CDRs H1
  abracd.pl purge.txt --purge --CDRs H1 H3 L1 --OUT test
  abracd.pl --CDRs ALL
  abracd.pl AB_loop_AbRaCD_H3.txt --raw --CDRs H3
  abracd.pl --raw --CDRs H1 H3
  abracd.pl --raw --CDRs ALL
  abracd.pl AbRaCD_input_H3.txt --model --CDRs H3

USAGE
}

sub fail {
    my ($message, $suggestion) = @_;
    log_event('ERROR', $message);
    log_event('SUGGESTION', $suggestion) if defined $suggestion;
    print STDERR "ERROR: $message\n\n";
    print STDERR "Suggested action: $suggestion\n\n" if defined $suggestion;
    usage(*STDERR);
    cleanup_temp_rcd_inputs();
    close $LOG_FH if $LOG_FH;
    exit 2;
}

sub write_unique_input_file {
    my ($directory, $stem, $kind, $lines_ref) = @_;
    for my $index (0 .. 9999) {
        my $suffix = $index ? "_${index}" : '';
        my $name = "${stem}_${kind}${suffix}.txt";
        my $path = File::Spec->catfile($directory, $name);
        my $flags = O_WRONLY | O_CREAT | O_EXCL;
        my $fh;
        if (!sysopen($fh, $path, $flags)) {
            next if $! == EEXIST;
            fail("Cannot create temporary RCD input '$path': $!",
                'Check write permissions for the input-list directory, then rerun: ' . $RERUN_COMMAND);
        }
        push @TEMP_RCD_INPUT_FILES, $path;
        print {$fh} @{$lines_ref} or fail("Could not write temporary RCD input '$path': $!",
            'Check available disk space and directory permissions, then rerun: ' . $RERUN_COMMAND);
        close $fh or fail("Could not close temporary RCD input '$path': $!",
            'Check directory permissions and rerun: ' . $RERUN_COMMAND);
        return File::Spec->rel2abs($path);
    }
    fail("Could not find an unused filename for the temporary H3 $kind input based on '$stem'.",
        'Remove or rename stale temporary files, then rerun: ' . $RERUN_COMMAND);
}

sub split_short_h3_input {
    my ($input_path) = @_;
    open my $input_fh, '<', $input_path or fail("Cannot read H3 RCD input '$input_path': $!",
        'Check that the H3 input list exists and is readable, then rerun: ' . $RERUN_COMMAND);
    my (@short_lines, @long_lines);
    my ($short_count, $long_count, $line_number) = (0, 0, 0);
    while (my $line = <$input_fh>) {
        ++$line_number;
        if ($line =~ /^\s*(?:#|$)/) {
            push @short_lines, $line;
            push @long_lines, $line;
            next;
        }
        my $parse_line = $line;
        $parse_line =~ s/^\s+//;
        my @fields = split /\s+/, $parse_line;
        fail("Malformed H3 RCD input at line $line_number in '$input_path'.",
            'Each entry must contain the PDB, start, end, chain, sequence, sequence length, and CDR fields.')
            unless @fields >= 7 && $fields[4] =~ /\S/ && $fields[5] =~ /^\d+$/ && uc($fields[6]) eq 'H3';
        my $sequence_length = length($fields[4]);
        fail("Sequence length mismatch at line $line_number in '$input_path' (sequence has $sequence_length residues, length field says $fields[5]).",
            'Correct the H3 RCD input list so its sequence and length fields agree.')
            unless $sequence_length == $fields[5];
        if ($sequence_length <= 7) {
            push @short_lines, $line;
            ++$short_count;
        } else {
            push @long_lines, $line;
            ++$long_count;
        }
    }
    close $input_fh or fail("Could not close H3 RCD input '$input_path': $!",
        'Check that the input list is readable, then rerun: ' . $RERUN_COMMAND);
    return unless $short_count;

    my $directory = dirname($input_path);
    my $stem = basename($input_path);
    $stem =~ s/\.txt\z//i;
    my $short_path = write_unique_input_file($directory, $stem, 'short', \@short_lines);
    my $long_path = write_unique_input_file($directory, $stem, 'long', \@long_lines);
    log_event('RCD INPUT', "H3 input split into $short_count short and $long_count long case(s).");
    return {
        short_path => $short_path,
        long_path => $long_path,
        short_count => $short_count,
        long_count => $long_count,
    };
}

sub run_rosetta_refinements {
    my ($cdr, $out, $rep, $rep_rcd, $rep_ros, $np, $model, $rcd_input_path, $refine_script) = @_;
    if ($rep_ros == 0) {
        print "Rosetta refinement skipped for $out because --repROS is 0.\n";
        log_event('ROSETTA', "Skipped refinement for $out because --repROS is 0.");
        return;
    }
    for my $ros (1 .. $rep_ros) {
        my $dir = $out;
        fail("RCD completed successfully but did not create directory '$dir'.",
            "Check RCD's output above and verify that the current directory is writable, then rerun: $RERUN_COMMAND")
            unless -d $dir;
        my $old = File::Spec->rel2abs('.');
        my $stage = "Rosetta refinement ($cdr, RCD repetition $rep/$rep_rcd, refinement $ros/$rep_ros)";
        $CURRENT_STAGE = $stage;
        chdir $dir or fail("Cannot enter output directory '$dir': $!",
            "Check directory permissions and rerun: $RERUN_COMMAND");
        my $suffix = "_${out}_ros${ros}";
        my $scenario = $model ? 'modeling' : 'native';
        my @refine_cmd = ($^X, $refine_script, $rcd_input_path, $scenario, $suffix, $np, '--energy ref2015');
        my $display_cmd = format_command(@refine_cmd);
        print "\n[$stage]\n> $display_cmd\n";
        log_event('ROSETTA', "Starting refinement $ros of $rep_ros in $dir; suffix: $suffix");
        log_event('COMMAND', $display_cmd);
        my $status = system { $refine_cmd[0] } @refine_cmd;
        chdir $old or fail("Cannot return to directory '$old': $!",
            "Restore access to the original working directory '$old' before starting another run.");
        if ($status != 0) {
            $CURRENT_STAGE = $stage;
            my $detail = $status == -1 ? "could not execute: $!"
                        : ($status & 127) ? "terminated by signal " . ($status & 127)
                        : "exit code " . ($status >> 8);
            fail("Stage '$stage' failed ($detail). Command: $display_cmd",
                "Review the Rosetta or PyRosetta error above, then rerun: $RERUN_COMMAND");
        }
        log_event('SUCCESS', "Rosetta refinement $ros of $rep_ros completed in $dir.");
        clean_model_outputs($dir) if $model;
    }
}

sub run_command {
    my ($stage, @cmd) = @_;
    my $options = ref($cmd[-1]) eq 'HASH' ? pop @cmd : {};
    my $quiet = $options->{quiet};
    $CURRENT_STAGE = $stage;
    my $display_cmd = format_command(@cmd);
    if ($quiet) {
        print "\n[$stage]\n";
    } else {
        print "\n[$stage]\n> $display_cmd\n";
    }
    log_event('COMMAND', $display_cmd);
    my $status = system { $cmd[0] } @cmd;
    if ($status == 0) {
        log_event('SUCCESS', 'Command completed successfully.');
        return 0;
    }
    if ($stage eq 'Chothia preprocessing' && $status != -1 && !($status & 127) && ($status >> 8) == 2) {
        log_event('PREPROCESS', 'Chothia preprocessing reported a critical CDR gap (exit code 2).');
        return $status;
    }
    my $detail = $status == -1 ? "could not execute: $!"
                : ($status & 127) ? "terminated by signal " . ($status & 127)
                : "exit code " . ($status >> 8);
    my $advice = $status == -1
        ? "Check that the executable and its dependencies are available, then rerun: $RERUN_COMMAND"
        : "Review the error reported by the program, correct its input or environment, then rerun: $RERUN_COMMAND";
    fail("Stage '$stage' failed ($detail). Command: $display_cmd", $advice);
}

sub take_option_value {
    my ($option, $args_ref, $default, $original_args_ref) = @_;
    fail("$option requires a value.",
        'Add a value. Corrected command: ' . corrected_option_command($option, $default, $original_args_ref))
        unless @$args_ref && $args_ref->[0] !~ /^--/;
    return shift @$args_ref;
}

sub corrected_option_command {
    my ($option, $value, $args_ref) = @_;
    my @corrected = @$args_ref;
    for my $index (0 .. $#corrected) {
        next unless $corrected[$index] eq $option;
        if ($index < $#corrected && $corrected[$index + 1] !~ /^--/) {
            $corrected[$index + 1] = $value;
        } else {
            splice @corrected, $index + 1, 0, $value;
        }
        return format_command($0, @corrected);
    }
    push @corrected, $option, $value;
    return format_command($0, @corrected);
}

sub command_without_option {
    my ($option, $args_ref) = @_;
    my @corrected;
    my $skip_value = 0;
    for my $arg (@$args_ref) {
        if ($skip_value) { $skip_value = 0; next; }
        if ($arg eq $option) {
            $skip_value = 1 if $option ne '--purge' && $option ne '--raw';
            next;
        }
        push @corrected, $arg;
    }
    return format_command($0, @corrected);
}

sub find_executable {
    my ($name) = @_;
    return $name if File::Spec->file_name_is_absolute($name) && -x $name;
    for my $dir (File::Spec->path()) {
        my $candidate = File::Spec->catfile($dir, $name);
        return $candidate if -f $candidate && -x $candidate;
    }
    return;
}

sub check_python_modules {
    my ($python, @modules) = @_;
    my $code = 'import importlib.util, sys; missing=[m for m in sys.argv[1:] if importlib.util.find_spec(m) is None]; print(" ".join(missing)); sys.exit(bool(missing))';
    my $status = system { $python } $python, '-c', $code, @modules;
    return $status == 0;
}

sub atomic_rewrite {
    my ($path, $lines_ref) = @_;
    my $tmp = "$path.tmp.$$";
    open my $out, '>', $tmp or die "Cannot write temporary file '$tmp': $!\n";
    print {$out} @$lines_ref;
    close $out or die "Cannot close temporary file '$tmp': $!\n";
    rename $tmp, $path or die "Cannot replace '$path' with '$tmp': $!\n";
}

sub trim_model_results {
    my ($path) = @_;
    fail("RCD results file not found: $path", 'Check that RCD completed successfully before model-output cleanup.')
        unless -f $path;
    open my $in, '<', $path or fail("Cannot read '$path': $!", 'Check file permissions and rerun the model workflow.');
    my (@kept, $found);
    while (my $line = <$in>) {
        if (!$found && $line =~ /^\s*RUN SPECIFICATIONS\s*$/) { $found = 1; }
        push @kept, $line if $found;
    }
    close $in;
    fail("Could not find RUN SPECIFICATIONS in $path", 'Keep the RCD results file for inspection; no lines were removed.')
        unless $found;
    atomic_rewrite($path, \@kept);
}

sub strip_rmsd_columns {
    my ($path) = @_;
    open my $in, '<', $path or fail("Cannot read Rosetta table '$path': $!", 'Check file permissions and rerun the model workflow.');
    my @output;
    while (my $line = <$in>) {
        $line =~ s/^\s+//;
        if ($line =~ /^$/ || $line =~ /^#[-]+\s*$/) { push @output, $line; next; }
        my @f = split /\s+/, $line;
        next unless @f;
        if ($f[0] =~ /^#loop$/i) {
            my @header = grep { defined } @f[0,1,2,4,8];
            @header = ('#Structure', '#loop', 'E_full', 'E_loop', 'E_KORP', 'E_full2')
                if basename($path) =~ /^rosetta/i;
            push @output, join("\t", @header) . "\n";
            next;
        }
        if (basename($path) =~ /_rosettaS?\.txt$/i && $f[0] eq '#' && defined($f[1]) && $f[1] =~ /^[-+]?\d/ ) {
            next;
        }
        my $has_structure = $f[0] !~ /^[-+]?\d+(?:\.\d+)?$/ && $f[0] ne '#'
            && $f[0] !~ /^#(?:Avg|Sig|Med)\(/;
        my @keep = $has_structure ? (0,1,2,3,5,9) : (0,1,2,4,8);
        if ($f[0] =~ /^#(?:Avg|Sig|Med)\(/) { @keep = (0,1,2,3,5,9); }
        my @row = map { $f[$_] } grep { $_ < @f } @keep;
        push @output, join("\t", @row) . "\n" if @row;
    }
    close $in;
    atomic_rewrite($path, \@output);
}

sub clean_model_outputs {
    my ($dir) = @_;
    my %files;
    $files{$_} = 1 for glob(File::Spec->catfile($dir, '*_rosetta.txt'));
    $files{$_} = 1 for glob(File::Spec->catfile($dir, '*_rosettaS.txt'));
    $files{$_} = 1 for glob(File::Spec->catfile($dir, 'rosetta*.txt'));
    fail("No Rosetta score tables found in '$dir'.", 'Check whether modeling refinement produced its score files.') unless %files;
    strip_rmsd_columns($_) for sort keys %files;
}

my @args = @ARGV;
my @original_args = @ARGV;
if (@args == 1 && ($args[0] eq '-h' || $args[0] eq '--help')) { usage(*STDOUT); exit 0; }

my $log_name_candidate;
for my $index (0 .. $#original_args - 1) {
    if ($original_args[$index] eq '--OUT'
        && $original_args[$index + 1] =~ /^[A-Za-z0-9._-]+$/) {
        $log_name_candidate = $original_args[$index + 1];
        last;
    }
}
$LOG_FILE = defined($log_name_candidate) ? "$log_name_candidate.log" : 'AbRaCD.log';
$RERUN_COMMAND = format_command($0, @original_args);
open $LOG_FH, '>', $LOG_FILE
    or fail("Cannot create log file '$LOG_FILE': $!", 'Check write permissions in the current working directory and rerun the command.');
$LOG_FH->autoflush(1);
log_event('START', "AbRaCD run started. Working directory: " . getcwd());
log_event('COMMAND', $RERUN_COMMAND);

my $input;
if (@args && $args[0] !~ /^--/) { $input = shift @args; }
my $raw = 0;
my $model = 0;
if (@args && ($args[0] eq '--raw' || $args[0] eq '--model')) {
    my $mode = shift @args;
    $raw = 1;
    $model = 1 if $mode eq '--model';
}
elsif (grep { $_ eq '--raw' || $_ eq '--model' } @args) {
    fail('--raw and --model must appear immediately after the input argument (or first when processing all PDB files).',
        'Use this order: abracd.pl input.txt --model --CDRs H3');
}
fail('--raw and --model cannot be used together.',
    'Keep --model only; it already selects the modeling input path.')
    if scalar(grep { $_ eq '--raw' } @original_args) && scalar(grep { $_ eq '--model' } @original_args);

my $purge = 0;
my ($out_name, $rep_rcd, $rep_ros, $nloops_override, $locobest_override) = (undef, 1, 1, undef, undef);
my @requested_cdrs;
my $np_requested = 0;
my %seen_option;
while (@args) {
    my $arg = shift @args;
    if ($arg eq '--purge') {
        fail('--purge may only be specified once.', 'Keep a single --purge option in the command.') if $seen_option{$arg}++;
        $purge = 1;
    }
    elsif ($arg eq '--CDRs') {
        fail('--CDRs may only be specified once.', 'Keep one --CDRs option and provide its CDR list once.') if $seen_option{$arg}++;
        while (@args && $args[0] !~ /^--/) {
            push @requested_cdrs, uc shift @args;
        }
        fail('--CDRs requires one or more CDR names, or ALL.',
            'Use --CDRs H1, --CDRs H1 H3, or --CDRs ALL.') unless @requested_cdrs;
    }
    elsif ($arg eq '--OUT') {
        fail('--OUT may only be specified once.', 'Keep one --OUT option and provide one output name.') if $seen_option{$arg}++;
        $out_name = take_option_value('--OUT', \@args, 'test', \@original_args);
    }
    elsif ($arg eq '--NP') {
        fail('--NP may only be specified once.', 'Keep one --NP option and provide one positive integer.') if $seen_option{$arg}++;
        $NP = take_option_value('--NP', \@args, '4', \@original_args);
        $np_requested = 1;
    }
    elsif ($arg eq '--NLOOPS') {
        fail('--NLOOPS may only be specified once.', 'Keep one --NLOOPS option and provide one positive integer.') if $seen_option{$arg}++;
        $nloops_override = take_option_value('--NLOOPS', \@args, '20000', \@original_args);
    }
    elsif ($arg eq '--LOCOBEST') {
        fail('--LOCOBEST may only be specified once.', 'Keep one --LOCOBEST option and provide one positive integer.') if $seen_option{$arg}++;
        $locobest_override = take_option_value('--LOCOBEST', \@args, '1000', \@original_args);
    }
    elsif ($arg eq '--repRCD') {
        fail('--repRCD may only be specified once.', 'Keep one --repRCD option and provide one positive integer.') if $seen_option{$arg}++;
        $rep_rcd = take_option_value('--repRCD', \@args, '2', \@original_args);
    }
    elsif ($arg eq '--repROS') {
        fail('--repROS may only be specified once.', 'Keep one --repROS option and provide a non-negative integer.') if $seen_option{$arg}++;
        $rep_ros = take_option_value('--repROS', \@args, '1', \@original_args);
    }
    else {
        fail("Unknown argument or invalid order: $arg",
            'Check the supported options in the usage above. Minimal valid command: abracd.pl --CDRs H3');
    }
}

log_event('INPUT', defined($input) ? "Input: $input" : 'Input: all PDB files in the current directory');
log_event('OPTIONS', "CDRs=" . (@requested_cdrs ? join(',', @requested_cdrs) : '(missing)')
    . "; raw=$raw; model=$model; purge=$purge; OUT=" . (defined($out_name) ? $out_name : '(default)')
    . "; repRCD=$rep_rcd; repROS=$rep_ros; NP=" . (defined($NP) ? $NP : '(missing)'));
print "Log file: $LOG_FILE\n";
fail('--CDRs is required.', 'Add --CDRs H3 or --CDRs ALL.') unless @requested_cdrs;
my @cdrs;
if (@requested_cdrs == 1 && $requested_cdrs[0] eq 'ALL') {
    @cdrs = qw(H1 H2 H3 L1 L2 L3);
} else {
    fail('ALL must be the only value after --CDRs.', 'Use --CDRs ALL by itself, or list specific CDRs such as --CDRs H1 H3.')
        if grep { $_ eq 'ALL' } @requested_cdrs;
    my %seen_cdr;
    for my $selected_cdr (@requested_cdrs) {
        fail("Invalid CDR: $selected_cdr. Accepted values: H1, H2, H3, L1, L2, L3.",
            'Use --CDRs followed by valid CDR names, for example --CDRs H1 H3.')
            unless exists $CDR_ARGS{$selected_cdr};
        fail("CDR '$selected_cdr' was listed more than once.", 'Remove duplicate CDR names from --CDRs.')
            if $seen_cdr{$selected_cdr}++;
        push @cdrs, $selected_cdr;
    }
}
fail('--OUT requires a non-empty name.',
    'Corrected command: ' . corrected_option_command('--OUT', 'test', \@original_args)) if defined($out_name) && !length($out_name);
fail('--OUT may contain only letters, numbers, dots, underscores, and hyphens.',
    'Corrected command: ' . corrected_option_command('--OUT', 'test', \@original_args))
    if defined($out_name) && $out_name !~ /^[A-Za-z0-9._-]+$/;
fail('--model requires exactly one CDR; multiple CDRs can overwrite design PDBs with identical names.',
    'Select one CDR, for example: abracd.pl input.txt --model --CDRs H3.')
    if $model && @cdrs != 1;
fail('--model requires one positional RCD input list.',
    'Pass the RCD list before --model, for example: abracd.pl input.txt --model --CDRs H3.')
    if $model && !defined($input);
fail('--model accepts only one RCD input .txt file.',
    'Pass an RCD-format .txt list as the positional input.')
    if $model && defined($input) && $input !~ /\.txt\z/i;
fail('--purge cannot be used with --raw or --model.',
    'Remove --purge, or use the crystallographic preprocessing workflow.') if $purge && $raw;
fail('--raw and --model are equivalent mode selectors; provide only one.',
    'Remove --raw or --model and keep the other mode selector.')
    if scalar(grep { $_ eq '--raw' } @original_args) && scalar(grep { $_ eq '--model' } @original_args);
for my $pair (['--repRCD', $rep_rcd]) {
    fail("$pair->[0] must be a positive integer.",
        'Use a positive integer, for example ' . corrected_option_command($pair->[0], '1', \@original_args))
        unless defined($pair->[1]) && $pair->[1] =~ /^\d+$/ && $pair->[1] > 0;
}
fail('--repROS must be a non-negative integer.',
    'Use 0 to skip Rosetta refinement, or a positive integer such as 1: '
        . corrected_option_command('--repROS', '1', \@original_args))
    unless defined($rep_ros) && $rep_ros =~ /^\d+$/;
for my $pair (['--NP', $NP], ['--NLOOPS', $nloops_override], ['--LOCOBEST', $locobest_override]) {
    next unless defined $pair->[1];
    my $default = $pair->[0] eq '--NP' ? '4'
                : $pair->[0] eq '--NLOOPS' ? '20000' : '1000';
    fail("$pair->[0] must be a positive integer.",
        'Use a positive integer, for example ' . corrected_option_command($pair->[0], $default, \@original_args))
        unless $pair->[1] =~ /^\d+$/ && $pair->[1] > 0;
}
$CURRENT_STAGE = 'input validation';
if ($purge && (!defined($input) || $input !~ /\.txt$/i)) {
    fail('With --purge, provide a .txt list with rows in this format: PDB H_chain L_chain.',
        "Provide a purge list as the input file, for example: abracd.pl purge.txt --purge --CDRs $cdrs[0]");
}
fail('--model and --raw cannot be used together.',
    'Keep --model only; it already selects the modeling input path.')
    if $model && scalar(grep { $_ eq '--raw' } @original_args);
if ($raw) {
    if (@cdrs == 1) {
        fail('With --raw and one selected CDR, provide its RCD input list as the positional input.',
            "Example: abracd.pl AbRaCD_input_$cdrs[0].txt --raw --CDRs $cdrs[0]") unless defined $input;
        fail("Raw RCD input must be a .txt file: $input",
            "Pass the RCD list for $cdrs[0], for example: abracd.pl AbRaCD_input_$cdrs[0].txt --raw --CDRs $cdrs[0]")
            unless $input =~ /\.txt$/i;
    } else {
        fail('With --raw and multiple CDRs or ALL, do not provide a positional input list.',
            'Remove the input filename. Example: abracd.pl --raw --CDRs H1 H3') if defined $input;
    }
} elsif (defined $input && $input !~ /\.(?:pdb|txt)$/i) {
    fail("Unsupported input file type: $input",
        'Use one PDB file or a .txt structure list. Example: abracd.pl structures.txt --CDRs ' . $cdrs[0]);
}
if (defined $input) {
    fail("Input file not found: $input", 'Check the input path and rerun: ' . $RERUN_COMMAND) unless -f $input;
}
if ($raw && @cdrs > 1) {
    my @missing_lists = map { "AbRaCD_input_${_}.txt" } grep { !-f "AbRaCD_input_${_}.txt" } @cdrs;
    fail('Missing raw RCD input list(s): ' . join(', ', @missing_lists),
        'Place every required AbRaCD_input_<CDR>.txt file in the current directory, or use a different --CDRs selection.')
        if @missing_lists;
}
$CURRENT_STAGE = 'preflight checks';
if (!$raw) {
    fail("Preprocessor not found: $PREPROCESS", 'Restore the AbRaCD Perl scripts and rerun: ' . $RERUN_COMMAND)
        unless -x $PREPROCESS || -f $PREPROCESS;
}
fail("Engineering preparation script not found: $ENGINEERING_PREP",
    'Restore ab_engineering_prep.pl beside abracd.pl and rerun: ' . $RERUN_COMMAND)
    if $model && !-f $ENGINEERING_PREP;
fail("Refinement script not found: $REFINE", 'Restore the AbRaCD Perl scripts and rerun: ' . $RERUN_COMMAND)
    if $rep_ros > 0 && !-f $REFINE;
for my $map ($CDR_MAP, $KORP_MAP, $DUN_MAP) {
    fail("Required RCD map not found: $map", 'Restore the required map file at the path shown above, then rerun: ' . $RERUN_COMMAND)
        unless -f $map;
}
$CURRENT_STAGE = 'preflight checks';
log_event('PREFLIGHT', 'Checking Python, maps, and RCD executables.');
my $python = find_executable('python3') || find_executable('python');
fail('Python 3 was not found on PATH; PyRosetta and mpi4py cannot be checked.',
    'Activate the Conda environment with AbRaCD dependencies, then rerun: ' . $RERUN_COMMAND) unless $python;
my @missing_modules;
push @missing_modules, 'PyRosetta (python module: pyrosetta)' unless check_python_modules($python, 'pyrosetta');
push @missing_modules, 'mpi4py (python module: mpi4py)' unless check_python_modules($python, 'mpi4py');
fail("Required Python dependencies are missing for '$python':\n  - " . join("\n  - ", @missing_modules),
    'Install the missing modules in the Python environment shown above, then rerun: ' . $RERUN_COMMAND)
    if @missing_modules;

my $rcd_executable;
my $use_mpi = 0;
my $mpi_launcher;
if ($np_requested) {
    $rcd_executable = find_executable('rcd_mpi_gcc')
        || find_executable('rcd_mpi_intel') || find_executable('rcd_mpi');
    fail('--NP was requested, but no MPI RCD executable was found on PATH (looked for rcd_mpi_gcc, rcd_mpi_intel, rcd_mpi).',
        'Add the RCD bin directory to PATH in the active environment, then rerun: ' . $RERUN_COMMAND)
        unless $rcd_executable;
    $use_mpi = 1;
    $mpi_launcher = find_executable('mpirun');
    fail('--NP was requested, but mpirun was not found on PATH.',
        'Activate or configure the MPI environment so mpirun is on PATH, then rerun: ' . $RERUN_COMMAND) unless $mpi_launcher;
} else {
    $rcd_executable = find_executable('rcd_gcc') || find_executable('rcd');
    unless ($rcd_executable) {
        $rcd_executable = find_executable('rcd_mpi_gcc')
            || find_executable('rcd_mpi_intel') || find_executable('rcd_mpi');
        fail('No RCD executable was found on PATH (looked for rcd_gcc, rcd, rcd_mpi_gcc, rcd_mpi_intel, rcd_mpi).',
            'Add the RCD bin directory to PATH in the active environment, then rerun: ' . $RERUN_COMMAND)
            unless $rcd_executable;
        $use_mpi = 1;
        $mpi_launcher = find_executable('mpirun');
        fail('Only an MPI RCD executable is available, but mpirun was not found on PATH.',
            'Activate or configure the MPI environment so mpirun is on PATH, then rerun: ' . $RERUN_COMMAND) unless $mpi_launcher;
    }
}
print "Preflight passed. Maps found; Python dependencies available; RCD executable: $rcd_executable\n";
print $use_mpi ? "RCD mode: MPI ($NP processes)\n" : "RCD mode: serial (default NP=1)\n";
log_event('PREFLIGHT', "Python=$python; RCD executable=$rcd_executable; MPI="
    . ($use_mpi ? "enabled ($NP processes), launcher=$mpi_launcher" : 'disabled'));
log_event('SETTINGS', "CDRs=" . join(',', @cdrs) . "; RCD repetitions=$rep_rcd; Rosetta repetitions per RCD=$rep_ros; energy=56; output prefix=" . (defined($out_name) ? $out_name : '(none)'));
log_event('PREFLIGHT', "Maps: CDR=$CDR_MAP; KORP=$KORP_MAP; Dunbrack=$DUN_MAP");
log_event('PREFLIGHT', 'All required programs, Python modules, and map files are available.');

my @pre_cmd = ($^X, $PREPROCESS);
push @pre_cmd, '--purge', $input if $purge;
push @pre_cmd, $input if defined($input) && !$purge;
if ($model) {
    my $stage = 'antibody engineering preparation';
    my @prep_cmd = ($^X, $ENGINEERING_PREP, $input);
    run_command($stage, @prep_cmd);
} elsif ($raw) {
    log_event('PREPROCESS', 'Skipped because --raw was specified.');
} else {
    my $pre_status = run_command('Chothia preprocessing', @pre_cmd);
    if ($pre_status && ($pre_status >> 8) == 2) {
        my @gap_cdrs = grep {
            my $gap_list = File::Spec->catfile(dirname(File::Spec->rel2abs($input // '.')), "AbRaCD_GAP_${_}.txt");
            -s $gap_list;
        } @cdrs;
        if (@gap_cdrs) {
            print "\nCritical structural gap detected for the selected CDR(s): ", join(', ', @gap_cdrs), ".\n";
            print "The engineering input list(s) were generated. Run the modeling workflow with:\n";
            for my $gap_cdr (@gap_cdrs) {
                my $gap_list = File::Spec->catfile(dirname(File::Spec->rel2abs($input // '.')), "AbRaCD_GAP_${gap_cdr}.txt");
                my @next = ($gap_list, '--model', '--CDRs', $gap_cdr);
                push @next, ('--NP', $NP) if $seen_option{'--NP'};
                push @next, ('--NLOOPS', $nloops_override) if $seen_option{'--NLOOPS'};
                push @next, ('--LOCOBEST', $locobest_override) if $seen_option{'--LOCOBEST'};
                push @next, ('--repRCD', $rep_rcd) if $seen_option{'--repRCD'};
                push @next, ('--repROS', $rep_ros) if $seen_option{'--repROS'};
                push @next, ('--OUT', $out_name) if defined $out_name;
                $next[0] = basename($gap_list);
                print "  ", join(' ', basename($0), @next), "\n";
            }
            log_event('ACTION REQUIRED', 'Critical gap detected; modeling command(s) printed for generated GAP list(s).');
            close $LOG_FH if $LOG_FH;
            exit 0;
        }
        log_event('WARNING', 'Preprocessing found a critical gap outside the selected CDR(s); continuing with selected input lists.');
    }
}

$CURRENT_STAGE = 'RCD input validation';
my %rcd_inputs;
for my $cdr (@cdrs) {
    my $list = $model
        ? File::Spec->catfile(dirname(File::Spec->rel2abs($input)), "AbRaCD_design_${cdr}.txt")
        : $raw && @cdrs == 1 ? $input : "AbRaCD_input_${cdr}.txt";
    $rcd_inputs{$cdr} = File::Spec->rel2abs($list);
}
my @missing_rcd_lists = grep { !-f $rcd_inputs{$_} } @cdrs;
if (@missing_rcd_lists) {
    my @missing_names = map { basename($rcd_inputs{$_}) } @missing_rcd_lists;
    my $advice = $raw && @cdrs == 1
        ? "Check that the supplied RCD list exists and contains valid entries for $cdrs[0], then rerun: $RERUN_COMMAND"
        : 'Run preprocessing first or place every required RCD list in the current directory: ' . join(', ', @missing_names);
    fail('Missing RCD input list(s): ' . join(', ', @missing_names), $advice);
}
log_event('RCD INPUT', 'Validated RCD input lists: ' . join(', ', map { basename($rcd_inputs{$_}) } @cdrs));

my @output_dirs;
for my $cdr (@cdrs) {
    my ($flag, $nokink, $dval, $default_nloops, $default_locobest) = @{$CDR_ARGS{$cdr}};
    my $nloops = defined($nloops_override) ? $nloops_override : $default_nloops;
    my $locobest = defined($locobest_override) ? $locobest_override : $default_locobest;
    my @cdr_extra_args = @{$CDR_EXTRA_ARGS{$cdr} || []};
    my $rcd_input_path = $rcd_inputs{$cdr};
    my $h3_split = $cdr eq 'H3' ? split_short_h3_input($rcd_input_path) : undef;
    my @deferred_rosetta;
    log_event('CDR START', "Starting CDR $cdr using " . basename($rcd_input_path)
        . "; flag=$flag; NLOOPS=$nloops; LOCOBEST=$locobest");

    for my $rep (1 .. $rep_rcd) {
        my $out = (defined($out_name) ? "${out_name}_" : '') . "${cdr}_rep${rep}";
        push @output_dirs, $out;
        log_event('RCD', "Starting RCD repetition $rep of $rep_rcd; output directory: $out");
        if ($h3_split) {
            my $short_results = File::Spec->catfile($out, 'results_shorts.txt');
            if ($h3_split->{short_count}) {
                fail("Refusing to overwrite existing short-loop results '$short_results'.",
                    'Choose a different --OUT prefix or move the existing result file before rerunning.') if -e $short_results;
                my @short_cmd = ($use_mpi ? ($mpi_launcher, '-np', $NP) : (), $rcd_executable,
                    $h3_split->{short_path}, '-r', '-e', 56, '--energy_file', $KORP_MAP,
                    '--nterm_rama', 0);
                push @short_cmd, '--bench' unless $model;
                push @short_cmd, ('-n', $nloops, '--loco_best', $locobest,
                    '-x', $DUN_MAP, '-d', 0.5, '--linear', '-t', 0.99, '-o', $out);
                run_command("RCD sampling ($cdr, repetition $rep/$rep_rcd)", @short_cmd, { quiet => 1 });
                trim_model_results(File::Spec->catfile($out, 'results.txt')) if $model;
                my $short_results_source = File::Spec->catfile($out, 'results.txt');
                fail("Short-loop RCD run did not create '$short_results_source'.",
                    "Check RCD's output above and verify the output directory is writable, then rerun: $RERUN_COMMAND")
                    unless -f $short_results_source;
                rename $short_results_source, $short_results
                    or fail("Could not preserve short-loop results as '$short_results': $!",
                        'Check directory permissions and available disk space, then rerun: ' . $RERUN_COMMAND);
            }
            if ($h3_split->{long_count}) {
                my @long_cmd = ($use_mpi ? ($mpi_launcher, '-np', $NP) : (), $rcd_executable,
                    $h3_split->{long_path}, '--ABmap', $CDR_MAP, $flag, @cdr_extra_args,
                    '--nokink_merge', $nokink, '--nterm_rama', 0, '-r', '-e', 56,
                    '--energy_file', $KORP_MAP, '-n', $nloops, '--loco_best', $locobest,
                    '-x', $DUN_MAP, '-d', $dval, '--linear', '-t', 0.99, '-o', $out);
                push @long_cmd, '--bench' unless $model;
                run_command("RCD sampling ($cdr, repetition $rep/$rep_rcd)", @long_cmd, { quiet => 1 });
                trim_model_results(File::Spec->catfile($out, 'results.txt')) if $model;
            }
        } else {
            my @rcd_cmd = ($use_mpi ? ($mpi_launcher, '-np', $NP) : (), $rcd_executable, $rcd_input_path,
                '--ABmap', $CDR_MAP, $flag, @cdr_extra_args, '--nokink_merge', $nokink,
                '--nterm_rama', 0, '-r', '-e', 56, '--energy_file', $KORP_MAP,
                '-n', $nloops, '--loco_best', $locobest, '-x', $DUN_MAP,
                '-d', $dval, '--linear', '-t', 0.99, '-o', $out);
            push @rcd_cmd, '--bench' unless $model;
            run_command("RCD sampling ($cdr, repetition $rep/$rep_rcd)", @rcd_cmd);
            trim_model_results(File::Spec->catfile($out, 'results.txt')) if $model;
        }

        if ($h3_split) {
            push @deferred_rosetta, [$out, $rep];
        } else {
            run_rosetta_refinements($cdr, $out, $rep, $rep_rcd, $rep_ros, $NP, $model, $rcd_input_path, $REFINE);
        }
    }
    if ($h3_split) {
        cleanup_temp_rcd_inputs();
        for my $deferred (@deferred_rosetta) {
            run_rosetta_refinements($cdr, $deferred->[0], $deferred->[1], $rep_rcd,
                $rep_ros, $NP, $model, $rcd_input_path, $REFINE);
        }
    }
}

print "\nWorkflow completed for CDRs " . join(', ', @cdrs) . ". RCD runs per CDR: $rep_rcd; Rosetta refinements per RCD run: $rep_ros.";
print " Rosetta refinement was skipped as requested (--repROS 0)." if $rep_ros == 0;
print "\n";
$CURRENT_STAGE = 'workflow complete';
log_event('COMPLETE', "Workflow completed. CDRs=" . join(',', @cdrs) . "; RCD runs per CDR=$rep_rcd; Rosetta refinements per RCD run=$rep_ros.");
log_event('ROSETTA', 'Rosetta refinement skipped by user (--repROS 0).') if $rep_ros == 0;
log_event('OUTPUT', "Output directories: " . join(', ', @output_dirs));
close $LOG_FH if $LOG_FH;
exit 0;
