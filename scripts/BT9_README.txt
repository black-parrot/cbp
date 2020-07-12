#*************************************************************
# (C) COPYRIGHT 2016 Samsung Electronics
#
#*************************************************************
* Dependency: Boost library
* BT9 File Format
  - Header section: 
    <key>: <value>
  - Node table
    NODE <id> <virtual_address> <physical_address> <opcode> <size> class: <br_class> behavior: <br_behavior> taken_cnt: <taken_cnt> tgt_cnt: <target_cnt> # mnemonic: "<assembly_code>"
  - Edge table
    EDGE <id> <src_id> <dest_id> <taken> <br_virt_target> <br_phy_target> <inst_cnt> traverse_cnt: <traverse_cnt>
  - Branch sequence list
    <edge_id>
    <edge_id>
    ...
    EOF
