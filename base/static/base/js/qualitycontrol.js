const f_sampleSelect = document.getElementById("f_sampleSelect");
const v_sampleSelect = document.getElementById("v_sampleSelect");
let sampletextnode;
let v_sampletextnode;
// proccess button click
async function processbtn(event) {
  event.preventDefault();
  toggleLoading(true, "processbutton");
  const csrftoken = getCookie("csrftoken");
  try {
    const process_result = await fetchAPI("/api/qualitycontrol", 0, csrftoken);
    if (process_result.error) {
      throw new Error(process_result.error.message);
    }
    // init modal
    initmodal();
    //   show btn
    const btnbox = document.getElementById("btnbox");
    btnbox.classList.remove("hidden");
    //   show image
    showimage(process_result.data);
    toggleLoading(false, "processbutton");
    // next page btn
    document.getElementById("nextbtn").classList.remove("hidden");
  } catch (error) {
    toggleLoading(false, "processbutton");
  }
}

function initmodal() {
  document.getElementById("filterform").reset();
  document.getElementById("previewimagebox").innerHTML = "";
  document.getElementById("confirmbtn").classList.add("hidden");
}
// filter method change
function showimage(result) {
  const imgbox = document.getElementById("imgbox");
  imgbox.textContent = "";
  const ul = document.createElement("ul");
  ul.classList.add(
    "space-y-1",
    "text-gray-500",
    "list-disc",
    "list-inside",
    "w-fit"
  );

  defaultsample = result["adata_results"][0].substring(
    0,
    result["adata_results"][0].indexOf(":")
  );
  sampletextnode = document.createTextNode(defaultsample);
  v_sampletextnode = document.createTextNode(defaultsample);
  const svalue = document.getElementById("sample");
  // const v_svalue = document.getElementById("v_sample");
  svalue.value = defaultsample;
  // v_svalue.value = defaultsample;
  f_sampleSelect.insertBefore(sampletextnode, f_sampleSelect.firstChild);
  // v_sampleSelect.insertBefore(v_sampletextnode, v_sampleSelect.firstChild);

  for (let i = 0; i < result.adata_results.length; i++) {
    const li = document.createElement("li");
    li.textContent = result["adata_results"][i];
    const img = document.createElement("img");
    img.src = `/media/qualitycontrol/origin/${
      result["save_image_names"][i]
    }?t=${Date.now()}`;
    img.classList.add("m-5", "size-3/4");
    li.appendChild(img);
    ul.appendChild(li);

    const sample_ul = document.getElementById("sampleul");
    const v_sample_ul = document.getElementById("v_sampleul");
    // selection add data
    var sample = result["adata_results"][i].substring(
      0,
      result["adata_results"][i].indexOf(":")
    );

    const sample_li = document.createElement("li");
    const v_sample_li = document.createElement("li");
    sample_li.innerHTML = `
              <a onclick="sampletextnodefix('${sample}')"
              class="block px-4 py-2 hover:bg-gray-100 dark:hover:bg-gray-600 dark:hover:text-white m-2"
          >${sample}</a>
      `;
    v_sample_li.innerHTML = `
            <a onclick="v_sampletextnodefix('${sample}')"
            class="block px-4 py-2 hover:bg-gray-100 dark:hover:bg-gray-600 dark:hover:text-white m-2"
        >${sample}</a>
    `;
    sample_ul.appendChild(sample_li);
    // v_sample_ul.appendChild(v_sample_li);
    console.log(sample_ul);
  }
  imgbox.appendChild(ul);
}
function sampletextnodefix(newText) {
  sampletextnode.nodeValue = newText;
  const sample = document.getElementById("sample");
  sample.value = newText;
}
function v_sampletextnodefix(newText) {
  v_sampletextnode.nodeValue = newText;
  const sample = document.getElementById("v_sample");
  sample.value = newText;
}
function methodtextnodefix(newText) {
  const methodtext = document.getElementById("methodtext");
  methodtext.textContent = newText;
  const fmethod = document.getElementById("fmethod");
  fmethod.value = newText;
  let upandlow = document.getElementById("upandlow");
  let lo_limit_label = document.getElementById("lo_limit_label");
  let up_limit_label = document.getElementById("up_limit_label");
  if (newText === "Manual") {
    upandlow.classList.remove("hidden");
    lo_limit_label.textContent = "Please enter lower limit:";
    up_limit_label.textContent = "Please enter upper limit:";
  } else if (newText === "Quantile") {
    upandlow.classList.remove("hidden");
    lo_limit_label.textContent = "Please enter quantile of lower limit:";
    up_limit_label.textContent = "Please enter quantile of upper limit:";
  } else {
    upandlow.classList.add("hidden");
  }
}

// preview
async function preview() {
  const previewimagebox = document.getElementById("previewimagebox");
  previewimagebox.textContent = "";
  const form = document.getElementById("filterform");
  const formData = new FormData(form);
  const csrftoken = getCookie("csrftoken");
  const preview_result = await fetchAPI("/api/preview", formData, csrftoken);
  // show preview image

  const ul = document.createElement("ul");
  ul.classList.add(
    "space-y-1",
    "text-gray-500",
    "list-disc",
    "list-inside",
    "w-fit"
  );
  const li = document.createElement("li");
  li.textContent = preview_result.data["adata_result"];
  const img = document.createElement("img");
  img.src = `/media/qualitycontrol/preview/${
    preview_result.data["save_image_names"]
  }?t=${Date.now()}`;
  img.classList.add("m-5", "size-3/4");
  li.appendChild(img);
  ul.appendChild(li);

  previewimagebox.appendChild(ul);
  //   show confirm btn
  document.getElementById("confirmbtn").classList.remove("hidden");
}
async function confirm() {
  const csrftoken = getCookie("csrftoken");
  const replaceimg_result = await fetchAPI("/api/replace", 0, csrftoken);

  showimage(replaceimg_result.data);
}
// before next page
async function presubmit(event) {
  const csrftoken = getCookie("csrftoken");
  event.preventDefault();
  confirm_result = await fetchAPI("/api/confirm", 0, csrftoken);
  window.location.href = event.target.href;
}
