use 5.006;
use strict;
use warnings FATAL => 'all';
use Module::Build;

# loading bin scripts #
opendir IN, "./bin/" or die $!;
my @binfiles = readdir IN;
my @binscripts = grep(/.pl$|.py$|.r$/i, @binfiles);
@binscripts = map{$_ = "./bin/" . $_} @binscripts;

my $builder = Module::Build->new(
    module_name         => 'ITEP_PopGen',
    license             => 'perl',
    dist_author         => q{Nick Youngblut <nyoungb2@illinois.edu>},
    dist_version_from   => 'lib/ITEP_PopGen.pm',
    release_status      => 'stable',
    configure_requires => {
        'Module::Build' => 0,
    },
    build_requires => {
        'Test::More' => 0,
    },
    requires => {
        #'ABC'              => 1.6,
        #'Foo::Bar::Module' => 5.0401,
    },
    add_to_cleanup     => [ 'ITEP_PopGen-*' ],
    create_makefile_pl => 'traditional',
    script_files => \@binscripts
);

$builder->create_build_script();

