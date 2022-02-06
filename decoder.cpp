#include "decoder.h"
#include <iostream>
#include <unordered_map>
#include <memory>
#include "huffman.h"
#include "fft.h"
#include <chrono>

class BitReader {
public:
    BitReader() = default;

    void SetInput(std::istream& input) {
        input_ = &input;
    }

    uint8_t GetByte() {
        curr_bit_ = -1;
        if (input_->eof()) {
            throw std::runtime_error("Unexpected end of input");
        }
        return input_->get();
    }

    uint16_t Peek() {
        return input_->peek();
    }

    uint16_t GetTwoBytes() {
        return (GetByte() << 8) | GetByte();
    }

    uint8_t GetFirstHalf(uint8_t byte) {
        return static_cast<uint8_t>(byte << 4) >> 4;
    }

    uint8_t GetSecondHalf(uint8_t byte) {
        return byte >> 4;
    }

    bool NextBit() {
        if (curr_bit_ == -1) {
            curr_byte_ = GetByte();
            if (curr_byte_ == 0xFF) {
                GetByte();
            }
            curr_bit_ = 7;
        }
        return curr_byte_ & (1 << curr_bit_--);
    }

    bool end = false;

    uint8_t curr_byte_ = 0, frame_ = 0;
    int curr_bit_ = -1;
    std::istream* input_;
};

static BitReader br;

enum BYTE {
    SECTION = 0xff,
    START = 0xd8,
    COMMENT = 0xfe,
    SPECIFIC1 = 0xe0,
    SPECIFIC2 = 0xe1,
    SPECIFIC3 = 0xe2,
    SPECIFIC4 = 0xec,
    SPECIFIC5 = 0xee,
    Q_TABLE = 0xdb,
    META_DATA = 0xc0,
    HUFFMAN = 0xc4,
    START_OF_SCANNING = 0xda,
    END = 0xd9
};

class Processor {
public:
    virtual ~Processor() = default;
    virtual void Process(Image&) = 0;
};

class ApplicationSpecific : public Processor {
public:
    void Process(Image&) override {
        uint16_t data_length = br.GetTwoBytes() - 2;
        for (uint16_t i = 0; i != data_length; ++i) {
            br.GetByte();
        }
    }
};

class Comment : public Processor {
public:
    void Process(Image& image) override {
        uint16_t data_length = br.GetTwoBytes() - 2;
        std::string comment;
        for (uint16_t i = 0; i != data_length; ++i) {
            comment += br.GetByte();
        }
        image.SetComment(comment);
    }
};

template <class T>
void SnakeFill(std::vector<T>& bytes, std::vector<std::vector<int>>& table) {
    size_t i = 0, j = 0, curr = 0;
    table[i][j] = bytes[curr++];
    while (curr != bytes.size()) {
        if (i == 0) {
            table[i][++j] = bytes[curr++];
            while (j != 0) {
                table[++i][--j] = bytes[curr++];
            }
        } else if (j == 0 && i != 7) {
            table[++i][j] = bytes[curr++];
            while (i != 0) {
                table[--i][++j] = bytes[curr++];
            }
        } else if (i == 7) {
            table[i][++j] = bytes[curr++];
            while (j != 7) {
                table[--i][++j] = bytes[curr++];
            }
        } else if (j == 7) {
            table[++i][j] = bytes[curr++];
            while (i != 7) {
                table[++i][--j] = bytes[curr++];
            }
        }
    }
}

class QuantizationTable : public Processor {
public:
    void Process(Image&) override {
        uint16_t data_length = br.GetTwoBytes() - 2;
        int cnt = 0;
        while (br.Peek() != BYTE::SECTION) {
            ++cnt;
            auto length_byte = br.GetByte();
            uint8_t values_length;
            if (br.GetSecondHalf(length_byte) == 0) {
                values_length = 1;
            } else if (br.GetSecondHalf(length_byte) == 1) {
                values_length = 2;
            } else {
                throw std::runtime_error("Wrong values length byte in DQT");
            }
            std::vector<int> bytes;
            bytes.reserve(64);
            for (size_t k = 0; k != 64; ++k) {
                int byte = values_length == 1 ? br.GetByte() : br.GetTwoBytes();
                bytes.push_back(byte);
            }
            q_tables[br.GetFirstHalf(length_byte)].resize(8, std::vector<int>(8));
            SnakeFill(bytes, q_tables[br.GetFirstHalf(length_byte)]);
        }
        if (data_length != 65 * cnt) {
            throw std::runtime_error("Bad DQT length");
        }
    }

    static std::unordered_map<uint8_t, std::vector<std::vector<int>>> q_tables;
};

std::unordered_map<uint8_t, std::vector<std::vector<int>>> QuantizationTable::q_tables = {};

struct Channel {
    uint8_t id;
    uint8_t h, v;
    uint8_t q_table_id;
    uint8_t cfh, cfv;
};

class BaselineMetaData : public Processor {
public:
    bool done = false;

    void Process(Image& image) override {
        if (done) {
            throw std::runtime_error("Second sof");
        }
        done = true;
        uint16_t data_length = br.GetTwoBytes();
        uint8_t precision = br.GetByte();
        uint16_t height = br.GetTwoBytes();
        uint16_t width = br.GetTwoBytes();
        uint8_t channels_num = br.GetByte();

        image.SetSize(width, height);

        uint8_t max_h = 1;
        uint8_t max_v = 1;
        for (uint8_t i = 0; i != channels_num; ++i) {
            uint8_t id = br.GetByte() - 1;
            uint8_t hv = br.GetByte();
            channels[id].id = id;
            channels[id].h = br.GetSecondHalf(hv);
            channels[id].v = br.GetFirstHalf(hv);
            max_h = std::max(max_h, channels[id].h);
            max_v = std::max(max_v, channels[id].v);
            channels[id].q_table_id = br.GetByte();
        }
        for (auto& ch : channels) {
            if (ch.second.h == 0 || ch.second.v == 0) {
                throw std::runtime_error("Baseline metadata");
            }
            ch.second.cfh = max_h / ch.second.h;
            ch.second.cfv = max_v / ch.second.v;
        }
    }

    static std::unordered_map<uint8_t, Channel> channels;
};

std::unordered_map<uint8_t, Channel> BaselineMetaData::channels = {};

struct HuffmanTable {
    HuffmanTree huffmanTree;
    uint8_t table_id;
};

class HuffmanBuild : public Processor {
public:
    void Process(Image&) override {
        uint16_t data_length = br.GetTwoBytes() - 3;
        int cnt = 2;
        while (br.Peek() != BYTE::SECTION) {
            uint8_t acdc_rock = br.GetByte();
            uint8_t type = br.GetSecondHalf(acdc_rock), table_id = br.GetFirstHalf(acdc_rock);
            std::vector<uint8_t> code_lengths, values;
            code_lengths.reserve(data_length);
            size_t values_num = 0;
            for (size_t i = 0; i != 16; ++i) {
                code_lengths.push_back(br.GetByte());
                values_num += code_lengths.back();
            }
            cnt += 17 + values_num;
            for (size_t i = 0; i != values_num; ++i) {
                values.push_back(br.GetByte());
            }
            if (type == 0) {
                dc_tables[table_id].table_id = table_id;
                dc_tables[table_id].huffmanTree.Build(code_lengths, values);
            } else if (type == 1) {
                ac_tables[table_id].table_id = table_id;
                ac_tables[table_id].huffmanTree.Build(code_lengths, values);
            } else {
                throw std::runtime_error("bad huffman table type");
            }
        }
        if (data_length + 3 != cnt) {
            throw std::runtime_error("huffman came wrong");
        }
    }

    static std::unordered_map<uint8_t, HuffmanTable> ac_tables, dc_tables;
};

class ScanStart : public Processor {
    struct Ids {
        uint8_t dc_huf_id;
        uint8_t ac_huf_id;
    };

public:
    void Process(Image& image) override {

        prev.resize(3);
        output.resize(64);
        DctCalculator calculator(8, &output, &output);
        uint16_t data_length = br.GetTwoBytes();
        uint8_t channels_num = br.GetByte();
        if (channels_num > 3) {
            throw std::runtime_error("Bad channels num");
        }
        std::vector<Ids> ids(channels_num);
        for (size_t i = 0; i != channels_num; ++i) {
            uint8_t id = br.GetByte();
            if (id < 1 || id > 3) {
                throw std::runtime_error("SOS bad id");
            }
            --id;
            uint8_t ac_dc_ids = br.GetByte();
            ids[id].ac_huf_id = br.GetFirstHalf(ac_dc_ids);
            ids[id].dc_huf_id = br.GetSecondHalf(ac_dc_ids);
            if (!HuffmanBuild::ac_tables.contains(ids[id].ac_huf_id) ||
                !HuffmanBuild::dc_tables.contains(ids[id].dc_huf_id)) {
                throw std::runtime_error("Bad huf table id");
            }
        }
        for (size_t i = 0; i != 3; ++i) {
            auto tmp = br.GetByte();
            if (i == 1 && tmp != 0x3f) {
                throw std::runtime_error("Bad progressive inf");
            }
        }
        std::vector<int> ch_amount(channels_num);
        for (size_t i = 0; i != channels_num; ++i) {
            ch_amount[i] = (BaselineMetaData::channels[i].cfh * BaselineMetaData::channels[i].cfv);
        }
        int y_mul = 8, x_mul = 8;
        if (channels_num == 1) {
            ch_amount[0] = 1;
        } else if (ch_amount[0] == 1 && ch_amount[1] == 4 && ch_amount[2] == 4) {
            ch_amount[0] = 4;
            ch_amount[1] = ch_amount[2] = 1;
            y_mul *= 2;
            x_mul *= 2;
        } else if (ch_amount[0] == 1 && ch_amount[1] == 2 && ch_amount[2] == 2) {
            ch_amount[0] = 2;
            ch_amount[1] = ch_amount[2] = 1;
            if (BaselineMetaData::channels[1].cfh == 1) {
                y_mul *= 2;
            } else if (BaselineMetaData::channels[1].cfh == 2) {
                x_mul *= 2;
            } else {
                throw std::runtime_error("baddddd");
            }
        }

        mcu_y.resize(ch_amount[0], std::vector<std::vector<int>>(8, std::vector<int>(8)));
        if (channels_num != 1) {
            mcu_cb.resize(ch_amount[1], std::vector<std::vector<int>>(8, std::vector<int>(8)));
            mcu_cr.resize(ch_amount[2], std::vector<std::vector<int>>(8, std::vector<int>(8)));
        }

        int blocks_in_row = image.Width() / x_mul + (image.Width() % x_mul != 0);
        int blocks_in_col = image.Height() / y_mul + (image.Height() % y_mul != 0);

        int whole = blocks_in_row * blocks_in_col;
        int y_shift = 0, x_shift = 0;
        int cnt = 0;

        BuildMapping();
        while (cnt != whole) {
            ReadPrep(mcu_y, ch_amount[0], 0, ids);
            if (channels_num != 1) {
                ReadPrep(mcu_cb, ch_amount[1], 1, ids);
                ReadPrep(mcu_cr, ch_amount[2], 2, ids);
            }
            QuantIDCTPrep(mcu_y, 0, calculator);
            if (channels_num != 1) {
                QuantIDCTPrep(mcu_cb, 1, calculator);
                QuantIDCTPrep(mcu_cr, 2, calculator);
            }
            SetRGBBlock(image, y_shift, x_shift, y_mul, x_mul);
            ++cnt;
            y_shift = (cnt / blocks_in_row) * y_mul;
            x_shift = (cnt % blocks_in_row) * x_mul;
        }
        br.end = true;
    }

    void ReadPrep(std::vector<std::vector<std::vector<int>>>& mcu_comp, int size, uint8_t id,
                  std::vector<Ids>& ids) {
        for (int i = 0; i != size; ++i) {
            ReadMatrix(id, ids, mcu_comp[i]);
            mcu_comp[i][0][0] += prev[id];
            prev[id] = mcu_comp[i][0][0];
        }
    }

    void QuantIDCTPrep(std::vector<std::vector<std::vector<int>>>& mcu_comp, uint8_t channel,
                       DctCalculator& calculator) {
        for (size_t i = 0; i != mcu_comp.size(); ++i) {
            QuantizationMul(channel, mcu_comp[i]);
            IDCT(mcu_comp[i], calculator);
        }
    }

    std::vector<double> output;

    void IDCT(std::vector<std::vector<int>>& mat, DctCalculator& calculator) {
        int k = 0;
        for (auto& row : mat) {
            for (auto el : row) {
                output[k++] = el;
            }
        }
        calculator.Inverse();
        k = 0;
        for (int i = 0; i != 8; ++i) {
            for (int j = 0; j != 8; ++j) {
                mat[i][j] = std::min(std::max(static_cast<int>(output[k++]) + 128, 0), 255);
            }
        }
    }

    std::vector<std::vector<std::vector<std::pair<int, std::pair<int, int>>>>> map;

    void BuildMapping() {
        map.resize(3, std::vector<std::vector<std::pair<int, std::pair<int, int>>>>(
                          16, std::vector<std::pair<int, std::pair<int, int>>>(
                                  16, std::pair<int, std::pair<int, int>>(-1, {0, 0}))));
        for (int id = 0; id != 3; ++id) {
            if (!BaselineMetaData::channels.contains(id)) {
                continue;
            }
            int h = BaselineMetaData::channels[id].cfh;
            int v = BaselineMetaData::channels[id].cfv;
            if (h < 1 || h > 2 || v < 1 || v > 2) {
                throw std::runtime_error("");
            }
            for (int i = 0; i != 16; ++i) {
                for (int j = 0; j != 16; ++j) {
                    int comp_num = 0;
                    if (mcu_y.size() == 4 && mcu_cb.size() == 1 && mcu_cr.size() == 1) {
                        comp_num = (i / 8) * 2 * (v & 1) + (j / 8) * (h & 1);
                    } else if (mcu_y.size() == 2 && mcu_cb.size() == 1 && mcu_cr.size() == 1) {
                        if (BaselineMetaData::channels[1].cfh == 1) {
                            comp_num = (i / 8) * (v & 1);
                        } else {
                            comp_num = (j / 8) * (h & 1);
                        }
                    }
                    map[id][i][j] = {comp_num,
                                     {(!(v & 1)) * (i / v) + (v & 1) * (i % 8),
                                      (!(h & 1)) * (j / h) + (h & 1) * (j % 8)}};
                }
            }
        }
    }

    void SetRGBBlock(Image& image, int y_shift, int x_shift, int y_size, int x_size) {
        for (int i = 0; i != y_size && i + y_shift < static_cast<int>(image.Height()); ++i) {
            for (int j = 0; j != x_size && j + x_shift < static_cast<int>(image.Width()); ++j) {
                SetPixel(image, i, j, y_shift, x_shift);
            }
        }
    }

    void SetPixel(Image& image, int ii, int jj, int y_shift, int x_shift) {
        auto [y_num, y_pos] = map.at(0).at(ii).at(jj);
        auto [cb_num, cb_pos] = map.at(1).at(ii).at(jj);
        auto [cr_num, cr_pos] = map.at(2).at(ii).at(jj);
        RGB pixel;
        double y_comp = mcu_y.at(y_num).at(y_pos.first).at(y_pos.second);
        double cb_comp =
            cb_num == -1 ? 128.0 : mcu_cb.at(cb_num).at(cb_pos.first).at(cb_pos.second);
        double cr_comp =
            cr_num == -1 ? 128.0 : mcu_cr.at(cr_num).at(cr_pos.first).at(cr_pos.second);

        pixel.r = int(y_comp + 1.402 * (cr_comp - 128.0));
        pixel.g = int(y_comp - 0.34414 * (cb_comp - 128.0) - 0.71414 * (cr_comp - 128.0));
        pixel.b = int(y_comp + 1.772 * (cb_comp - 128.0));
        pixel.r = std::min(std::max(0, pixel.r), 255);
        pixel.g = std::min(std::max(0, pixel.g), 255);
        pixel.b = std::min(std::max(0, pixel.b), 255);
        image.SetPixel(ii + y_shift, jj + x_shift, pixel);
    }

    void QuantizationMul(uint8_t channel, std::vector<std::vector<int>>& mat) {
        if (!QuantizationTable::q_tables.contains(BaselineMetaData::channels[channel].q_table_id)) {
            throw std::runtime_error("Qtable mismatch");
        }
        auto& table = QuantizationTable::q_tables[BaselineMetaData::channels[channel].q_table_id];
        if (table.size() != mat.size()) {
            throw std::runtime_error("abacaba");
        }
        for (size_t i = 0; i != table.size(); ++i) {
            for (size_t j = 0; j != table.size(); ++j) {
                mat[i][j] *= table[i][j];
            }
        }
    }

    uint8_t DiveToTerminate(HuffmanTree& ht) {
        int ans;
        while (!ht.Move(br.NextBit(), ans)) {
        }
        return ans;
    }

    int GetNumByLen(uint8_t len) {
        bool inverted = false;
        int ans = 0;
        for (size_t i = 0; i != len; ++i) {
            auto bit = br.NextBit();
            if (i == 0 && !bit) {
                inverted = true;
            }
            ans |= bit << (len - 1 - i);
        }
        if (inverted) {
            ans -= (1 << len) - 1;
        }
        return ans;
    }

    void ReadMatrix(uint8_t id, std::vector<Ids>& ids, std::vector<std::vector<int>>& mcu_comp) {
        HuffmanTree& dc_tree = HuffmanBuild::dc_tables[ids[id].dc_huf_id].huffmanTree;
        HuffmanTree& ac_tree = HuffmanBuild::ac_tables[ids[id].ac_huf_id].huffmanTree;

        uint8_t dc_len = DiveToTerminate(dc_tree);
        int dc_coef = GetNumByLen(dc_len);
        std::vector<double> nums(64);
        nums[0] = dc_coef;
        uint8_t ac_byte = DiveToTerminate(ac_tree);
        int cnt = 1;
        while (ac_byte != 0) {
            uint8_t zeros = br.GetSecondHalf(ac_byte);
            cnt += zeros;
            uint8_t ac_len = br.GetFirstHalf(ac_byte);
            int ac_coef = GetNumByLen(ac_len);
            nums[cnt++] = ac_coef;
            if (cnt >= 64) {
                break;
            }
            ac_byte = DiveToTerminate(ac_tree);
        }
        SnakeFill(nums, mcu_comp);
    }

    std::vector<std::vector<std::vector<int>>> mcu_y, mcu_cb, mcu_cr;
    std::vector<int> prev;
};

std::unordered_map<uint8_t, HuffmanTable> HuffmanBuild::ac_tables = {};
std::unordered_map<uint8_t, HuffmanTable> HuffmanBuild::dc_tables = {};

void ReadSectionByte() {
    if (br.GetByte() != BYTE::SECTION) {
        throw std::runtime_error("Section expected");
    }
}

class Cleaner {
public:
    ~Cleaner() {
        br = {};
        HuffmanBuild::dc_tables.clear();
        HuffmanBuild::ac_tables.clear();
        BaselineMetaData::channels.clear();
        QuantizationTable::q_tables.clear();
    }
};

Image Decode(std::istream& input) {
    Cleaner cleaner;
    const std::unordered_map<uint8_t, std::shared_ptr<Processor>> processors{
        {BYTE::COMMENT, std::make_shared<Comment>()},
        {BYTE::SPECIFIC1, std::make_shared<ApplicationSpecific>()},
        {BYTE::SPECIFIC2, std::make_shared<ApplicationSpecific>()},
        {BYTE::SPECIFIC3, std::make_shared<ApplicationSpecific>()},
        {BYTE::SPECIFIC4, std::make_shared<ApplicationSpecific>()},
        {BYTE::SPECIFIC5, std::make_shared<ApplicationSpecific>()},
        {BYTE::Q_TABLE, std::make_shared<QuantizationTable>()},
        {BYTE::META_DATA, std::make_shared<BaselineMetaData>()},
        {BYTE::HUFFMAN, std::make_shared<HuffmanBuild>()},
        {BYTE::START_OF_SCANNING, std::make_shared<ScanStart>()}};

    Image ans;
    br.SetInput(input);
    ReadSectionByte();
    if (br.GetByte() != BYTE::START) {
        throw std::runtime_error("Bad image start");
    }
    while (!br.end) {
        ReadSectionByte();
        processors.at(br.GetByte())->Process(ans);
    }
    ReadSectionByte();
    if (br.GetByte() != BYTE::END) {
        throw std::runtime_error("Bad end");
    }
    return ans;
}
