# These example shows that the VCF output line can be edited. (Thanks to Shane McCarthy)

{
    tag      => 'FORMAT/GQ',
    name     => 'MinSampleGQ',
    desc     => 'Genotypes set to . for samples with GQ < 10', 
    apply_to => 'all',
    test     => sub {
        my $i = 8;
        for my $gq (@$MATCH)
        {
            $i++;
            next unless ($gq<10); 
            my @format = split(/:/,$$RECORD[$i]);
            $format[0] = $format[0] =~ /\// ? "./." : ".";
            $$RECORD[$i] = join(":",@format);
        }
        return $PASS;
    },
},
