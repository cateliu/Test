transcript on
if {[file exists rtl_work]} {
	vdel -lib rtl_work -all
}
vlib rtl_work
vmap work rtl_work

vlog -vlog01compat -work work +incdir+D:/Altera/Documents/Key_Cnt/rtl {D:/Altera/Documents/Key_Cnt/rtl/counter.v}
vlog -vlog01compat -work work +incdir+D:/Altera/Documents/Key_Cnt/rtl {D:/Altera/Documents/Key_Cnt/rtl/uart_byte_tx.v}
vlog -vlog01compat -work work +incdir+D:/Altera/Documents/Key_Cnt/rtl {D:/Altera/Documents/Key_Cnt/rtl/key_time_cnt.v}
vlog -vlog01compat -work work +incdir+D:/Altera/Documents/Key_Cnt/rtl {D:/Altera/Documents/Key_Cnt/rtl/ctrl.v}

vlog -vlog01compat -work work +incdir+D:/Altera/Documents/Key_Cnt/prj/testbench {D:/Altera/Documents/Key_Cnt/prj/testbench/key_time_cnt_tb.v}

vsim -t 1ps -L altera_ver -L lpm_ver -L sgate_ver -L altera_mf_ver -L altera_lnsim_ver -L cycloneive_ver -L rtl_work -L work -voptargs="+acc"  key_time_cnt_tb

add wave *
view structure
view signals
run -all
