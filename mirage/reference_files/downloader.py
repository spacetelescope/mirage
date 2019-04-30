#! /usr/bin /env python

"""Download the reference files needed to run Mirage. Extract and unzip
files, and place into the proper directory. Inform the user how to set
the MIRAGE_DATA encironment variable.
"""
import os
import requests
import shutil
import tarfile
import gzip

from mirage.utils.utils import ensure_dir_exists


NIRCAM_REFFILES_URL = [('https://stsci.box.com/shared/static/6eomezd68n3surgqut8if6gy8l6lf3xk.gz', 'nircam_reference_files.tar.gz')]
NIRISS_REFFILES_URL = [('https://stsci.box.com/shared/static/evlv7vxszgmiff3zdmdxnol6u6h8pa9j.gz', 'niriss_reference_files.tar.gz')]
FGS_REFFILES_URL = [('https://stsci.box.com/shared/static/ia5z21m69tb08hd5zpv01c43g0px3gfm.gz', 'fgs_reference_files.tar.gz')]

NIRCAM_CR_LIBRARY_URL = [('https://stsci.box.com/shared/static/4cw7wmsqw9qhdozl4owz0tmr6ozusfqr.gz', 'mirage_nircam_cr_library.tar.gz')]
NIRISS_CR_LIBRARY_URL = [('https://stsci.box.com/shared/static/uxyb08cjkf1i7yd4fhryrhi6dr4da9pg.gz', 'mirage_niriss_cr_library.tar.gz')]
FGS_CR_LIBRARY_URL = [('https://stsci.box.com/shared/static/d5oswszqbwt6i027g6ue3usi47dmyign.gz', 'mirage_fgs_cr_library.tar.gz')]



NIRCAM_INDIVIDUAL_PSF_URLS = [('https://stsci.box.com/shared/static/4rw5p9dd8ofa13pnmgfgz1svr8qustcv.gz', 'nircam_a1_subpix_grid_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/ngx0mfo6pppe9zbrvzh28n90fe6v4d6i.gz', 'nircam_a2_subpix_grid_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/nc2tt2w4cbylwy0ny7bxme37xtnd26cm.gz', 'nircam_a3_subpix_grid_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/27bkgppw3ajyhfvl6nnv2wt9mlcn0ae9.gz', 'nircam_a4_subpix_grid_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/gg48crkex2afjqb11replcg0gs2bm5cv.gz', 'nircam_a5_subpix_grid_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/tg2pjk829oc73ijl9kda0j4ccagg75ce.gz', 'nircam_b1_subpix_grid_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/7vmp33fb6hwcyesl39z7v7rcgpq4ls5x.gz', 'nircam_b2_subpix_grid_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/55yf2n72l0xlc6u6t6igv5e8u8zcnunq.gz', 'nircam_b3_subpix_grid_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/5gbse7nhpwnfptpip5bci3uyt0sate4y.gz', 'nircam_b4_subpix_grid_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/b3fepnceqznfmwct9sszoc2lzxl8iaku.gz', 'nircam_b5_subpix_grid_webbpsf_library.tar.gz')]
NIRISS_INDIVIDUAL_PSF_URLS = [('https://stsci.box.com/shared/static/29z7ounfn61gjwi3kiu04hc0f8qgaot9.gz', 'niriss_subpix_grid_f090w_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/i5yhe4zbzxawnah8csdpp4qtk9dmao1t.gz', 'niriss_subpix_grid_f115w_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/s1j1kpx0aes1ntf4qrlyl1mhitj0sz2u.gz', 'niriss_subpix_grid_f140m_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/cc5kjqp78uvtrmlpd0jk6llobajkhavr.gz', 'niriss_subpix_grid_f150w_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/g0lwik43s7erbpsm6r5slq1739bnjob9.gz', 'niriss_subpix_grid_f158m_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/70utbxqfmi474hxb34uawjnrujilj6hb.gz', 'niriss_subpix_grid_f200w_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/jorfqdq56u73wlnj9nk0r34o36ruc8ky.gz', 'niriss_subpix_grid_f277w_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/edafbfzxc4iijr5rzfkziasu64quh6lb.gz', 'niriss_subpix_grid_f356w_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/oqsriubvc0knwtrsoavxj8gi7psz2yck.gz', 'niriss_subpix_grid_f380m_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/i6c28jii139x8dnzxy6bp3fnm25npnvw.gz', 'niriss_subpix_grid_f430m_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/fjep4n3dhw2uavt0uk07dfk65dww1zp8.gz', 'niriss_subpix_grid_f444w_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/3cg3397siwiirv64hpez5uc9oyaey5d8.gz', 'niriss_subpix_grid_f480m_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/rk38t0psx4kmqx5wxqsmlrntzusm4aod.gz', 'niriss_nrm_subpix_grid_f277w_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/x2c99ivjzze0ixdywtoczoyi438mw76m.gz', 'niriss_nrm_subpix_grid_f380m_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/lwwei3atsbo64iz10c3blru5fnbof69u.gz', 'niriss_nrm_subpix_grid_f430m_webbpsf_library.tar.gz'),
                              ('https://stsci.box.com/shared/static/7mxdonrlx9qz1mxmfgjzfg6bqwc5fzi3.gz', 'niriss_nrm_subpix_grid_f480m_webbpsf_library.tar.gz')]

FGS_INDIVIDUAL_PSF_URLS = [('https://stsci.box.com/shared/static/3g8f3i0w24l4yqu0bpin5uei7e5or9uh.gz', 'fgs_subpix_grid_webbpsf_library.tar.gz')]

NIRCAM_GRIDDED_PSF_URLS = []
NIRISS_GRIDDED_PSF_URLS = []
FGS_GRIDDED_PSF_URLS = []

NIRCAM_RAW_DARK_URLS = [('https://stsci.box.com/shared/static/pctjgthruh86ctr6ww9bgzccjlzc0xt3.gz', 'NRCNRCA1-DARK-60082202011_1_481_SE_2016-01-09T00h03m58_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/df2bxy3iwenot0ykti06xqwnn1j9k0ii.gz', 'NRCNRCA1-DARK-60090213141_1_481_SE_2016-01-09T02h53m12_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/rv8x4snxwqznymy92h40c70f7c7k8zh7.gz', 'NRCNRCA1-DARK-60090604481_1_481_SE_2016-01-09T06h52m47_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/v0h1ndnnbgiar4wcx1ynl0ueli4ksz57.gz', 'NRCNRCA1-DARK-60091005411_1_481_SE_2016-01-09T10h56m36_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/c1u5iwte7brjlm2dhqfibbw6ylqyshjv.gz', 'NRCNRCA1-DARK-60091434481_1_481_SE_2016-01-09T15h50m45_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/60mpcp81s7nzmjmb6cbnjvznxqpqbja8.gz', 'NRCNRCA2-DARK-60082224241_1_482_SE_2016-01-09T00h10m36_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/9zbro9it89h6ddipmfjxsakf3mrleb7x.gz', 'NRCNRCA2-DARK-60090235001_1_482_SE_2016-01-09T04h17m03_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/fe96k5h9polxacaewwit01n2ufutndpq.gz', 'NRCNRCA2-DARK-60090635511_1_482_SE_2016-01-09T07h05m19_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/3uw9i5q2uwgsol6q5qb7uzf4c3l295ro.gz', 'NRCNRCA2-DARK-60091030561_1_482_SE_2016-01-09T11h03m17_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/deuxfqva8vd32oblmz3d0wu8683doigs.gz', 'NRCNRCA2-DARK-60091457131_1_482_SE_2016-01-09T15h50m45_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/8mx9ziefvjil6ptsr52bqatg69ssr832.gz', 'NRCNRCA3-DARK-60082245481_1_483_SE_2016-01-09T00h04m26_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/e6pm1zg59k7fdx36x9cxl0kalvd5ewwo.gz', 'NRCNRCA3-DARK-60090321241_1_483_SE_2016-01-09T04h17m10_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/qf90unkzkfb7gv0fhg3n1kkf2hwfc6l8.gz', 'NRCNRCA3-DARK-60090656591_1_483_SE_2016-01-09T07h31m27_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/8p348edl8dh5wr7g2sqnj2tlxae3us9o.gz', 'NRCNRCA3-DARK-60091052561_1_483_SE_2016-01-09T11h28m06_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/zoq314y0jzanpyd0am0cu3f9phoklmzv.gz', 'NRCNRCA3-DARK-60091522581_1_483_SE_2016-01-09T16h30m34_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/sv5npcnyvkfshdsjgi0ooen5edwxsg2e.gz', 'NRCNRCA4-DARK-60082307391_1_484_SE_2016-01-09T00h04m08_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/zgjvynxemydgiu435tpjh7kt30ywk9wp.gz', 'NRCNRCA4-DARK-60090259591_1_484_SE_2016-01-09T04h16m50_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/q8s7f2icgdtnzjd4a2qat9gmm30vqfhj.gz', 'NRCNRCA4-DARK-60090720041_1_484_SE_2016-01-09T07h58m26_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/ya8o0881reg9jncshug5b91mug3trxlk.gz', 'NRCNRCA4-DARK-60091117441_1_484_SE_2016-01-09T11h52m23_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/zifn5ogf7glksnddhptbqbvxioxtw2dm.gz', 'NRCNRCA4-DARK-60091548131_1_484_SE_2016-01-09T16h30m50_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/tjkyvudkgwxsehw6xuizz8bhduy2dkwu.gz', 'NRCNRCALONG-DARK-60082329041_1_485_SE_2016-01-09T00h04m16_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/9kp7zwv88baoqb1v2ajx2x96yubkv1ml.gz', 'NRCNRCALONG-DARK-60090344021_1_485_SE_2016-01-09T04h16m42_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/14osssajxnyfhwnhkmk60ob8u25ymiuq.gz', 'NRCNRCALONG-DARK-60090746381_1_485_SE_2016-01-09T08h21m48_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/f2ntx77wg82padmj4u9wjuzeoxdiy4mq.gz', 'NRCNRCALONG-DARK-60091140151_1_485_SE_2016-01-09T14h23m49_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/ghn59ko04b2msldt089iffc1v3ie7agi.gz', 'NRCNRCALONG-DARK-60091611271_1_485_SE_2016-01-09T17h16m35_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/lfyonc2kuvevakzr10yrgcgjk7rcm6o8.gz', 'NRCNRCB1-DARK-60082356471_1_486_SE_2016-01-09T02h47m00_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/jhhbkxlw47aey5ut93o5xbnk8rrrrcdj.gz', 'NRCNRCB1-DARK-60090405201_1_486_SE_2016-01-09T05h33m56_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/padrxurzkvg4dda9ombf7fto2poqcwtl.gz', 'NRCNRCB1-DARK-60090807581_1_486_SE_2016-01-09T08h48m11_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/bhzkf4uex4yqa2sq0pzt1srhtq6liv13.gz', 'NRCNRCB1-DARK-60091205311_1_486_SE_2016-01-09T14h30m08_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/04m3kxsbtehf7ndvxut70jz814g6lr77.gz', 'NRCNRCB1-DARK-60091636021_1_486_SE_2016-01-09T17h16m13_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/wsifqidhgwixkpjpmqy9i4anhagf361h.gz', 'NRCNRCB2-DARK-60090021181_1_487_SE_2016-01-09T02h51m54_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/dubivztvsafigzm7oghyi7l7hvns3az5.gz', 'NRCNRCB2-DARK-60090427541_1_487_SE_2016-01-09T05h33m14_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/70nmsc1op9mq1uybceh3m5zw547otjgl.gz', 'NRCNRCB2-DARK-60090830131_1_487_SE_2016-01-09T08h59m50_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/qo3b6icrzn3ztd2t27nq1acexb321fqx.gz', 'NRCNRCB2-DARK-60091230011_1_487_SE_2016-01-09T14h23m47_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/jlr9f0uiajyoaf5j3bqas4iqgk7mgxla.gz', 'NRCNRCB2-DARK-60091735131_1_487_SE_2016-01-09T18h09m45_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/b0bfeqkaeqgyblzhl39ufa9t0p8b4vh1.gz', 'NRCNRCB3-DARK-60090043151_1_488_SE_2016-01-09T02h53m21_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/nx7c1ugce6factspliruya3z3t7chkmx.gz', 'NRCNRCB3-DARK-60090451471_1_488_SE_2016-01-09T05h33m25_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/k8uuqkouz5465t5m5238rf327ehly7ol.gz', 'NRCNRCB3-DARK-60090852451_1_488_SE_2016-01-09T09h35m03_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/4wd55isjm8s2zdmqsf5lc9fdiwzopwfh.gz', 'NRCNRCB3-DARK-60091254111_1_488_SE_2016-01-09T14h23m58_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/tr722w87zdvyee35sou62ipjkhgd6z2b.gz', 'NRCNRCB3-DARK-60091757401_1_488_SE_2016-01-09T18h40m55_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/runlc3olxq3b6uw56cf8gq0ck5mx227p.gz', 'NRCNRCB4-DARK-60090118261_1_489_SE_2016-01-09T02h46m53_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/5v9xqtswzuvprgshp2zpunbp5h32izkt.gz', 'NRCNRCB4-DARK-60090513431_1_489_SE_2016-01-09T05h57m50_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/r7hf20hxg2mamfuleaneq8zz138hvwog.gz', 'NRCNRCB4-DARK-60090914351_1_489_SE_2016-01-09T09h52m02_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/ked1azkg05ia2qledmpm29mkiktinn3y.gz', 'NRCNRCB4-DARK-60091316411_1_489_SE_2016-01-09T14h23m38_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/ninym1sr7yrgb9tjmfczkb1vh6v2an0t.gz', 'NRCNRCB4-DARK-60091822061_1_489_SE_2016-01-09T18h53m02_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/y9psja8ol4kqvkqjx4z5hlo7cbqrpqp7.gz', 'NRCNRCBLONG-DARK-60090141241_1_490_SE_2016-01-09T02h46m50_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/7jj24l6i3k5yy23gclqvvs1gy8cqeiix.gz', 'NRCNRCBLONG-DARK-60090535381_1_490_SE_2016-01-09T06h17m51_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/56h3wdjro3ici53tbt6nivjc0ntxx8f3.gz', 'NRCNRCBLONG-DARK-60090939281_1_490_SE_2016-01-09T10h22m25_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/8wq3o6n8bgh448zpgz7801lnndzh24da.gz', 'NRCNRCBLONG-DARK-60091338491_1_490_SE_2016-01-09T15h46m43_level1b_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/bks6n8q83ocfremqn99jt8c5kn8niyid.gz', 'NRCNRCBLONG-DARK-60101408431_1_490_SE_2016-01-10T15h01m09_level1b_uncal.fits.gz')]
NIRCAM_LINEARIZED_DARK_URLS = [('https://stsci.box.com/shared/static/dchmjjmepngpdo4xrabk8d5h1q63m3hf.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA1-DARK-60082202011_1_481_SE_2016-01-09T00h03m58_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/z5fp0yhy7hydb3m59yevgxif8b9rrg8f.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA1-DARK-60090213141_1_481_SE_2016-01-09T02h53m12_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/mqbwnwx8vjibgd50hahtd4i5qvc1mo0p.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA1-DARK-60090604481_1_481_SE_2016-01-09T06h52m47_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/dc0y3k95jo5pd3yhskig4bimvgfjkopl.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA1-DARK-60091005411_1_481_SE_2016-01-09T10h56m36_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/8j4jbzis927j74m6eizy90n0l5ukbd2g.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA1-DARK-60091434481_1_481_SE_2016-01-09T15h50m45_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/6h3uau1yr2sdb16xyn7dbfv83gmq1iaf.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA2-DARK-60082224241_1_482_SE_2016-01-09T00h10m36_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/re2lt2feyvyypj8pjn47tug8qlk1hzfm.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA2-DARK-60090235001_1_482_SE_2016-01-09T04h17m03_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/fszoeeviszj91li1v9amcfce68g7f5hz.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA2-DARK-60090635511_1_482_SE_2016-01-09T07h05m19_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/3gukyemqcjcqnt39f7s2sxyjlzcr7mk3.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA2-DARK-60091030561_1_482_SE_2016-01-09T11h03m17_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/via9puv3co1yeqg6vrqpsftxheer6m6d.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA2-DARK-60091457131_1_482_SE_2016-01-09T15h50m45_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/cnjlguana504waizgw15ez19lfgx3mwo.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA3-DARK-60082245481_1_483_SE_2016-01-09T00h04m26_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/2epggrnq3misb3zn8e7iys9sxxhgzrqt.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA3-DARK-60090321241_1_483_SE_2016-01-09T04h17m10_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/3dpeo1kwn4orhpgxdybm3qo830rlcvm9.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA3-DARK-60090656591_1_483_SE_2016-01-09T07h31m27_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/7nin41qwumemceb6ma3xurqycg4dde1w.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA3-DARK-60091052561_1_483_SE_2016-01-09T11h28m06_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/6yt605vznl90sjo41ny8lgwtdyhxzpu7.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA3-DARK-60091522581_1_483_SE_2016-01-09T16h30m34_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/5abqwlhau5cvgqtf7kinpgstqwx5tf7d.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA4-DARK-60082307391_1_484_SE_2016-01-09T00h04m08_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/1h86n0fcnmex8cfrjsl56hzdcjkbqaox.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA4-DARK-60090259591_1_484_SE_2016-01-09T04h16m50_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/wttboeuk25eaxovslkm1psckitmtdeuh.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA4-DARK-60090720041_1_484_SE_2016-01-09T07h58m26_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/acspa28skwaagjje74sairln2k61s37l.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA4-DARK-60091117441_1_484_SE_2016-01-09T11h52m23_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/l9171f1e8x3g8ichtu3w6gze6g8oo7ik.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCA4-DARK-60091548131_1_484_SE_2016-01-09T16h30m50_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/8vrtbwiquazo3v0gis2ibd677ro5jw75.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCALONG-DARK-60082329041_1_485_SE_2016-01-09T00h04m16_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/mk15yiw25jc6ucsg8f8hl9xya55wvcqz.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCALONG-DARK-60090344021_1_485_SE_2016-01-09T04h16m42_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/f8ikd8ursxxmizh5i88rzz1frvg7n002.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCALONG-DARK-60090746381_1_485_SE_2016-01-09T08h21m48_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/7zzdnsgrg5kddp067hv6umpztbzido1s.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCALONG-DARK-60091140151_1_485_SE_2016-01-09T14h23m49_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/2djmrufvnbu014i1uilyphqbogcplj03.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCALONG-DARK-60091611271_1_485_SE_2016-01-09T17h16m35_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/bebwkwderw0y4qhg342z5joyvbm3ddmu.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB1-DARK-60082356471_1_486_SE_2016-01-09T02h47m00_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/kbicf1i4lp4lv951lakxn7i5mak40jg2.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB1-DARK-60090405201_1_486_SE_2016-01-09T05h33m56_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/f7egrgol987bkev77cy5c63ljy4r0hv5.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB1-DARK-60090807581_1_486_SE_2016-01-09T08h48m11_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/ba0lgcegd8h3xormuc41wle2i4fr1p1e.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB1-DARK-60091205311_1_486_SE_2016-01-09T14h30m08_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/44chr1i2gyojtke34pboj7nyfh2mtmne.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB1-DARK-60091636021_1_486_SE_2016-01-09T17h16m13_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/8f01ir2s6lt6itisoypw3t2owvyajwno.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB2-DARK-60090021181_1_487_SE_2016-01-09T02h51m54_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/isrohgs3ekpt9a2qubojqwgt2vkefyxa.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB2-DARK-60090427541_1_487_SE_2016-01-09T05h33m14_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/vh5pv2as2nfzh3rgutuwdr9kw4iz0kwt.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB2-DARK-60090830131_1_487_SE_2016-01-09T08h59m50_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/p6ibgezka3vkvr5bsjv52xlvius925ff.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB2-DARK-60091230011_1_487_SE_2016-01-09T14h23m47_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/n94tmd71cgabs3dswsdodj1p65uwz0cw.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB2-DARK-60091735131_1_487_SE_2016-01-09T18h09m45_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/sb1xugk7cdq6zhqkb65ifbl5y1h3l8ob.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB3-DARK-60090043151_1_488_SE_2016-01-09T02h53m21_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/m5zqf0fqq5ycw9n0v5vageukk1v4wl9j.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB3-DARK-60090451471_1_488_SE_2016-01-09T05h33m25_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/gqkru1s8zwvrps3l4t31d1prk6f5r2b9.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB3-DARK-60090852451_1_488_SE_2016-01-09T09h35m03_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/upk0cqi9aqivj300ui2s8non2a5u3jlb.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB3-DARK-60091254111_1_488_SE_2016-01-09T14h23m58_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/lh7595afpjqvriwaajrqd9sm3i54zxhl.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB3-DARK-60091757401_1_488_SE_2016-01-09T18h40m55_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/pmmr6ozzxo4zgbxv7oqwg1xbyrze67mu.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB4-DARK-60090118261_1_489_SE_2016-01-09T02h46m53_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/nrui4vu3jgptopw812mgjf0pxjnjh95k.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB4-DARK-60090513431_1_489_SE_2016-01-09T05h57m50_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/faxxltcjpm45g2fugvd8ewu4kodm8a3z.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB4-DARK-60090914351_1_489_SE_2016-01-09T09h52m02_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/ll9n77ca9i8dz7ae41b9ae4hc1jlh5x6.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB4-DARK-60091316411_1_489_SE_2016-01-09T14h23m38_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/59u1h4sxthx2brm6ukaacbbcpa36wuhm.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCB4-DARK-60091822061_1_489_SE_2016-01-09T18h53m02_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/mhy8m4sdhcfv4t5bkqbtn0sr8lzplsxo.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCBLONG-DARK-60090141241_1_490_SE_2016-01-09T02h46m50_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/1t03er9udcj0fjacx1huhyx5qj8chjzv.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCBLONG-DARK-60090535381_1_490_SE_2016-01-09T06h17m51_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/qmssc39jl1wm3bv2biiawha67xc17eiv.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCBLONG-DARK-60090939281_1_490_SE_2016-01-09T10h22m25_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/q0sbu5lekvl1gk5u9wrvkhm3e74m7vfq.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCBLONG-DARK-60091338491_1_490_SE_2016-01-09T15h46m43_uncal.fits.gz'),
                               ('https://stsci.box.com/shared/static/0hf9nopvq0vk16oselokbnqo80wmf788.gz', 'Linearized_Dark_and_SBRefpix_NRCNRCBLONG-DARK-60101408431_1_490_SE_2016-01-10T15h01m09_uncal.fits.gz')]

NIRISS_RAW_DARK_URLS = [('https://stsci.box.com/shared/static/1pr9cfwx2d8r6iju9afmhsylowtwq3x4.gz', 'NISNIRISSDARK-153451235_11_496_SE_2015-12-11T16h05m20_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/cbnz7he6l9nva1lhrmn8pvz42hfb2ke5.gz', 'NISNIRISSDARK-153451235_12_496_SE_2015-12-11T16h23m51_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/xxwjvy4w6hvmllkkzhdpydwk4hy9p997.gz', 'NISNIRISSDARK-153451235_13_496_SE_2015-12-11T16h42m52_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/0a8xl3e6w45rw5s5xlnk4m0nkhx2nrgd.gz', 'NISNIRISSDARK-153451235_14_496_SE_2015-12-11T17h01m50_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/jjax31gw1bxp4iiaeen4gk869yic3a9k.gz', 'NISNIRISSDARK-153451235_15_496_SE_2015-12-11T17h20m40_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/rw04slw0ib53giuey7igysm6ge2hbzet.gz', 'NISNIRISSDARK-153451235_16_496_SE_2015-12-11T17h40m30_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/n3geovm7mxn31f4piziqb89htu4ak1bu.gz', 'NISNIRISSDARK-153451235_17_496_SE_2015-12-11T17h59m52_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/8ddnqcphujkncpg6niy1wgq5d5otlth1.gz', 'NISNIRISSDARK-153451235_18_496_SE_2015-12-11T18h16m31_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/umituldjs916c5yhl0at8si29mg10qmq.gz', 'NISNIRISSDARK-153451235_19_496_SE_2015-12-11T18h36m32_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/mbhyo8eqescy4v77gwznqd0aajk72wx2.gz', 'NISNIRISSDARK-153451235_20_496_SE_2015-12-11T18h53m52_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/kh8ep1ecazy9fmw4we36kxg51e75hhzi.gz', 'NISNIRISSDARK-172500017_13_496_SE_2017-09-07T04h48m22_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/1dc1q3f942gzjr76r81837tt6uyfh5qk.gz', 'NISNIRISSDARK-172500017_14_496_SE_2017-09-07T05h06m42_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/dg2niytahewqpsepto1tl23q7u03b230.gz', 'NISNIRISSDARK-172500017_15_496_SE_2017-09-07T05h28m22_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/99preol28qzhg7ctt0r4vsf8nrblcitt.gz', 'NISNIRISSDARK-172500017_16_496_SE_2017-09-07T05h47m42_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/agtzclaccw93gobr76kqygvbs3xth0n2.gz', 'NISNIRISSDARK-172500017_17_496_SE_2017-09-07T06h09m02_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/dr4hkbw15397iy93wmqcqr1mcxijdss5.gz', 'NISNIRISSDARK-172500017_18_496_SE_2017-09-07T06h29m12_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/ebo2m6kaaqwvxu5rtnpiky6amt6k8733.gz', 'NISNIRISSDARK-172500017_19_496_SE_2017-09-07T06h49m52_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/y0db7ckqj4y7cfhmkegphi3arxacews0.gz', 'NISNIRISSDARK-172500017_20_496_SE_2017-09-07T07h09m22_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/3t3by35z9jhz1uz7su85rxwkw4j3812m.gz', 'NISNIRISSDARK-172500017_21_496_SE_2017-09-07T07h29m52_dms_uncal.fits.gz'),
                        ('https://stsci.box.com/shared/static/5udxxbpnkda7fjpgc52wd9lv125gwt0i.gz', 'NISNIRISSDARK-172500017_22_496_SE_2017-09-07T07h50m32_dms_uncal.fits.gz')]
NIRISS_LINEARIZED_DARK_URLS = [('https://stsci.box.com/shared/static/0lb7f57m9iakhiuc6p22kc3okgtcj10q.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_17_496_SE_2015-12-11T17h59m52_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/cduow17sa0cejtdm28rh48um9mll9co0.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_16_496_SE_2015-12-11T17h40m30_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/gvx1mjhqpmbtvpln3aiojfigsfflb2v3.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_15_496_SE_2015-12-11T17h20m40_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/ymee55dyp1dnarofrgk9inlgoo5c9lko.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_14_496_SE_2015-12-11T17h01m50_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/by38gaqd2w8x3l0fmuox76dc26wj3ry2.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_13_496_SE_2015-12-11T16h42m52_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/vif9ha3z70mf4x4ln70464o093qalbjs.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_12_496_SE_2015-12-11T16h23m51_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/nte0e1k6mtmym50kqe5wlngg6ambm6j2.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_11_496_SE_2015-12-11T16h05m20_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/tzlrfwwp7cwvnstc7h8h4zvtfidpqeyh.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_18_496_SE_2015-12-11T18h16m31_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/96pvv9yo3liew5tduij4t9egydwh7m81.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_19_496_SE_2015-12-11T18h36m32_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/1b8v1zugso34uu6gpvplcq19qclhazta.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_20_496_SE_2015-12-11T18h53m52_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/j76aa9m6qob2n4iix728o50ug738m0r3.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_13_496_SE_2017-09-07T04h48m22_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/0ncv03i5uhp6e97y8ofb31j3flzd5e87.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_14_496_SE_2017-09-07T05h06m42_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/w1lcdxo927qr6eq0pqx5eiit8iz5c3ae.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_15_496_SE_2017-09-07T05h28m22_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/3xf32m9k9t9daukblqpcjje2yj9mx10d.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_16_496_SE_2017-09-07T05h47m42_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/clud1dkupuxkkf77c96qbhdyfn0ukeop.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_17_496_SE_2017-09-07T06h09m02_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/nhrpau4kuxs0pm93orng1k6n6mh29m06.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_18_496_SE_2017-09-07T06h29m12_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/vjel9s4x44nqn55d8dsk9zauwq96tz3h.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_19_496_SE_2017-09-07T06h49m52_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/o52djwdwksvnjkl7xt0ndlnjpfa6tlms.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_20_496_SE_2017-09-07T07h09m22_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/857olzwa1e51b05uz2bx6c13hzsufilz.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_21_496_SE_2017-09-07T07h29m52_dms_uncal_linear_dark_prep_object.fits.gz'),
                               ('https://stsci.box.com/shared/static/isnevg8m2j58qrb72u19pp43k5rpksy1.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_22_496_SE_2017-09-07T07h50m32_dms_uncal_linear_dark_prep_object.fits.gz')]

# For testing
print('URL list shortened for testing')
NIRISS_RAW_DARK_URLS = [('https://stsci.box.com/shared/static/1pr9cfwx2d8r6iju9afmhsylowtwq3x4.gz', 'NISNIRISSDARK-153451235_11_496_SE_2015-12-11T16h05m20_dms_uncal.fits.gz')]
NIRISS_LINEARIZED_DARK_URLS = [('https://stsci.box.com/shared/static/0lb7f57m9iakhiuc6p22kc3okgtcj10q.gz', 'Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_17_496_SE_2015-12-11T17h59m52_dms_uncal_linear_dark_prep_object.fits.gz')]


FGS_RAW_DARK_URLS = [('https://stsci.box.com/shared/static/2nhm4pajg1d3b3vmj8p5wtsevxq41qsj.gz', '29722_1x88_FGSF03512-D-NR-G2-5339214947_1_498_SE_2015-12-05T22h27m19_dms_uncal.fits.gz'),
                     ('https://stsci.box.com/shared/static/yq0pvyur651h8v1fz9t5wronvtbhv76b.gz', '29782_1x88_FGSF03872-PAR-5340074326_1_498_SE_2015-12-06T12h22m47_dms_uncal.fits.gz'),
                     ('https://stsci.box.com/shared/static/72byn3psj6g3oawh3kp6fx5szfy1ahfs.gz', '29813_1x88_FGSF037221-MR-2-5340161743_1_498_SE_2015-12-06T16h45m10_dms_uncal.fits.gz'),
                     ('https://stsci.box.com/shared/static/kqp2rmgff8esq2dyi5zeduaayunyhgtz.gz', '30632_1x88_FGSF03511-D-NR-G1-5346180117_1_497_SE_2015-12-12T19h00m12_dms_uncal.fits.gz'),
                     ('https://stsci.box.com/shared/static/j3nlgllcmsmi8wtxo5kewzb7wvrefq0i.gz', '30670_1x88_FGSF03511-D-NR-G2-5346181816_1_498_SE_2015-12-12T21h31m01_dms_uncal.fits.gz'),
                     ('https://stsci.box.com/shared/static/of8p9yk3jeo8hmi3k91wi4gk9kol8vlo.gz', '30742_1x88_FGSF03871-PAR-5347035139_1_498_SE_2015-12-13T05h23m30_dms_uncal.fits.gz'),
                     ('https://stsci.box.com/shared/static/54ggqqbf2ltnfy0efed5d2othk4sf0rv.gz', '30749_1x88_FGSF03881-PAR-5347043800_1_497_SE_2015-12-13T09h02m01_dms_uncal.fits.gz'),
                     ('https://stsci.box.com/shared/static/2iilllsjo54yj7g8m6d931u19aforcbf.gz', '30829_1x88_FGSF037111-G1NRNC-5347151640_1_497_SE_2015-12-13T16h28m38_dms_uncal.fits.gz')]
FGS_LINEARIZED_DARK_URLS = [('https://stsci.box.com/shared/static/6y0rsqgongmyp9cffqwd6k2k3ivl6pjf.gz', '29722_1x88_FGSF03512-D-NR-G2-5339214947_1_498_SE_2015-12-05T22h27m19_dms_uncal_linearized.fits.gz'),
                            ('https://stsci.box.com/shared/static/t7xyd85wcvvxgzaeh26wlmstmp437g39.gz', '29782_1x88_FGSF03872-PAR-5340074326_1_498_SE_2015-12-06T12h22m47_dms_uncal_linearized.fits.gz'),
                            ('https://stsci.box.com/shared/static/y5dr624wwc3rf7iadc51dblw9vn5bavf.gz', '29813_1x88_FGSF037221-MR-2-5340161743_1_498_SE_2015-12-06T16h45m10_dms_uncal_linearized.fits.gz'),
                            ('https://stsci.box.com/shared/static/xg35tec3ohaeihmp1qpiktur2y9lnfem.gz', '30632_1x88_FGSF03511-D-NR-G1-5346180117_1_497_SE_2015-12-12T19h00m12_dms_uncal_linearized.fits.gz'),
                            ('https://stsci.box.com/shared/static/0rkd1vzphgqzjl8yinusck4c57dm5nov.gz', '30670_1x88_FGSF03511-D-NR-G2-5346181816_1_498_SE_2015-12-12T21h31m01_dms_uncal_linearized.fits.gz'),
                            ('https://stsci.box.com/shared/static/4xe1mmap18612b57joctwo20t38ald95.gz', '30742_1x88_FGSF03871-PAR-5347035139_1_498_SE_2015-12-13T05h23m30_dms_uncal_linearized.fits.gz'),
                            ('https://stsci.box.com/shared/static/jq206o2pm5ydud1xidxyurcf74zwuzb2.gz', '30749_1x88_FGSF03881-PAR-5347043800_1_497_SE_2015-12-13T09h02m01_dms_uncal_linearized.fits.gz'),
                            ('https://stsci.box.com/shared/static/kai7ms09pvbfm5zdr079gss27ki3j6qd.gz', '30829_1x88_FGSF037111-G1NRNC-5347151640_1_497_SE_2015-12-13T16h28m38_dms_uncal_linearized.fits.gz')]


def download_file(url, file_name, output_directory='./'):
    """Download into the current working directory the
    file from Box given the direct URL

    Parameters
    ----------
    url : str
        URL to the file to be downloaded

    Returns
    -------
    download_filename : str
        Name of the downloaded file
    """
    download_filename = os.path.join(output_directory, file_name)
    if not os.path.isfile(download_filename):
        print('Downloading: {}'.format(file_name))
        with requests.get(url, stream=True) as response:
            if response.status_code != 200:
                raise RuntimeError("Wrong URL - {}".format(url))
            with open(download_filename, 'wb') as f:
                for chunk in response.iter_content(chunk_size=2048):
                    if chunk:
                        f.write(chunk)
    return download_filename


def download_reffiles(directory, instrument='all', psf_version='subpix', dark_type='linearized'):
    """Download tarred and gzipped reference files. Expand, unzip and
    organize into the necessary directory structure such that Mirage
    can use them.

    Parameters
    ----------
    directory : str
        Directory into which the reference files are placed. This will
        be the directory set to the MIRAGE_DATA environment variable

    instrument : str
        If ``all``: download all files. If the name of an individual
        instrument, download only the data for that instrument

    psf_version : str
        If ``gridded``: download new PSf library
        If ``something else``: download old PSF library

    dark_type : str
        Type of dark current files to download. Options are:
        'linearized': download linearize dark current ramps
        'raw': download raw dark current ramps
        'both': download both raw and linearized dark current ramps
    """
    file_list = get_file_list(instrument.lower(), psf_version.lower(), dark_type.lower())

    # Download everything first
    for file_info in file_list:
        file_url, filename = file_info
        download_file(file_url, filename, directory)
        local_file = os.path.join(directory, filename)

    # Now untar/organize. This way if the download is interrupted, it can
    # pixk up where it left off, since no downloaded files will have been
    # moved yet.
    for file_info in file_list:
        file_url, filename = file_info
        local_file = os.path.join(directory, filename)
        # Unzip and untar file
        if 'tar.gz' in local_file:
            print('Unzipping/extracting {}'.format(filename))
            file_object = tarfile.open(name=local_file, mode='r:gz')
            file_object.extractall(path=directory)
        else:
            # Darks need to be unzipped into the correct directory
            if 'linearized' in filename.lower():
                cal = 'linearized'
            else:
                cal = 'raw'

            # Determine directory
            if 'NRCNRC' in filename:
                det_str = filename.split('NRCNRC')[1].split('-')[0]
                if 'LONG' in det_str:
                    det_str.replace('LONG', '5')
                darks_dir = os.path.join(directory, 'mirage_data', 'nircam', 'darks')
                sub_directory = os.path.join(darks_dir, cal, det_str)
            elif 'niriss' in filename.lower():
                darks_dir = os.path.join(directory, 'mirage_data', 'niriss', 'darks')
                sub_directory = os.path.join(darks_dir, cal)
                #sub_directory = os.path.join(directory, 'mirage_data', 'niriss', 'darks', cal)
            elif 'fgs' in filename:
                darks_dir = os.path.join(directory, 'mirage_data', 'fgs', 'darks')
                sub_directory = os.path.join(darks_dir, cal)

            # Create the directory if it does not yet exist
            ensure_dir_exists(darks_dir)
            ensure_dir_exists(os.path.join(darks_dir, cal))
            ensure_dir_exists(sub_directory)

            final_location = os.path.join(sub_directory, filename)

            # Move the zipped file into the correct subdirectory
            if not os.path.isfile(final_location):
                print('Moving {} to {}'.format(filename, sub_directory))
                shutil.move(local_file, final_location)

            # Unzip
            unzipped_filename = final_location.replace('.gz', '')
            if not os.path.isfile(unzipped_filename):
                print('Unzipping {}'.format(filename))
                unzip_file(final_location)

    print(('Mirage reference files downloaded and extracted. Before '
           'using Mirage, be sure to set the MIRAGE_DATA environment '
           'variable to point to {}/mirage_data'.format(directory)))
    print('\n In bash: ')
    print('export MIRAGE_DATA="{}"'.format(os.path.join(directory, 'mirage_data')))


def get_file_list(instruments, library_version, dark_current):
    """Collect the list of URLs corresponding to the Mirage reference
    files to be downloaded

    Parameters
    ----------
    instruments : list
        List of instrument names for which to download data

    library_version : str
        Version of the PSF librarry to download.
        If ``gridded``: download PSf library that uses griddedPSFModels
        If ``something else``: download PSF library composed of
        individual fits files

    dark_current : str
        Type of dark current files to download. Options are:
        'linearized': download linearize dark current ramps
        'raw': download raw dark current ramps
        'both': download both raw and linearized dark current ramps

    Returns
    -------
    urls : list
        List of tuples. Each tuple contains:
        (URL for downloading, name of file)
    """
    urls = []
    instrument_names = [name.strip().lower() for name in instruments.split(',')]

    if 'all' in instrument_names:
        instrument_names = ['nircam', 'niriss', 'fgs']

    for instrument_name in instrument_names:
        # NIRCam
        if instrument_name.lower() == 'nircam':
            urls.extend(NIRCAM_REFFILES_URL)
            urls.extend(NIRCAM_CR_LIBRARY_URL)

            if library_version == 'gridded':
                urls.extend(NIRCAM_GRIDDED_PSF_URLS)
            else:
                urls.extend(NIRCAM_INDIVIDUAL_PSF_URLS)

            if dark_current in ['linearized', 'both']:
                urls.extend(NIRCAM_LINEARIZED_DARK_URLS)
            elif dark_current in ['raw', 'both']:
                urls.extend(NIRCAM_RAW_DARK_URLS)

        # NIRISS
        elif instrument_name.lower() == 'niriss':
            urls.extend(NIRISS_REFFILES_URL)
            urls.extend(NIRISS_CR_LIBRARY_URL)

            if library_version == 'gridded':
                urls.extend(NIRISS_GRIDDED_PSF_URLS)
            else:
                urls.extend(NIRISS_INDIVIDUAL_PSF_URLS)

            if dark_current in ['linearized', 'both']:
                urls.extend(NIRISS_LINEARIZED_DARK_URLS)
            elif dark_current in ['raw', 'both']:
                urls.extend(NIRISS_RAW_DARK_URLS)

        # FGS
        elif instrument_name.lower() == 'fgs':
            urls.extend(NIRISS_REFFILES_URL)
            urls.extend(NIRISS_CR_LIBRARY_URL)

            if library_version == 'gridded':
                urls.extend(FGS_GRIDDED_PSF_URLS)
            else:
                urls.extend(FGS_INDIVIDUAL_PSF_URLS)

            if dark_current in ['linearized', 'both']:
                urls.extend(FGS_LINEARIZED_DARK_URLS)
            elif dark_current in ['raw', 'both']:
                urls.extend(FGS_RAW_DARK_URLS)
    return urls


def unzip_file(filename):
    """Unzip a file

    Parameters
    ----------
    filename : str
        Name of file to unzip

    dir_name : str
        Directory into which the file is unzipped
    """
    if not os.path.isfile(filename):
        print('File {} does not exist.'.format(filename))
    with gzip.open(filename, 'rb') as file_in:
        unzipped_name = filename.replace('.gz', '')
        with open(unzipped_name, 'wb') as file_out:
            shutil.copyfileobj(file_in, file_out)
