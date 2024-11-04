// grab cookie
function getCookie(name) {
  let cookieValue = null;
  if (document.cookie && document.cookie !== "") {
    const cookies = document.cookie.split(";");
    for (let i = 0; i < cookies.length; i++) {
      const cookie = cookies[i].trim();
      if (cookie.substring(0, name.length + 1) === name + "=") {
        cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
        break;
      }
    }
  }
  return cookieValue;
}
// fetchapi
async function fetchAPI(url, formData, csrftoken) {
  try {
    const response = await fetch(url, {
      method: "POST",
      body: formData,
      headers: {
        "X-CSRFToken": csrftoken,
      },
    });

    const data = await response.json();

    if (!response.ok) {
      throw { message: data, status: response.status };
    }

    return { data, status: response.status };
  } catch (error) {
    return { error: error.message || "未知錯誤", status: error.status || 500 };
  }
}
// grab steps
async function grabsteps(url) {
  try {
    const response = await fetch(url);
    if (!response.ok) {
      return null;
    }
    const data = await response.json();
    return data;
  } catch (error) {
    return null;
  }
}
// del steps
async function deleteSteps() {
  let step = document.getElementById("steps").value;
  let steps = {
    createadata: 1,
    qualitycontrol: 2,
    normalizationmerge: 5,
    pca: 7,
    clustering: 8,
    identifythegates: 14,
    phenotyping: 16,
    spatialanalysis: 22,
    neighborhood: 23,
  };
  let url = `/delSteps/${steps[step]}/`;
  try {
    const response = await fetch(url);
    if (!response.ok) {
      return null;
    }
    window.location.reload();
  } catch (error) {
    return null;
  }
}
function downloadImages(folder) {
  const csrftoken = getCookie("csrftoken");
  fetch("/download_image/", {
    method: "POST",
    headers: {
      "Content-Type": "application/x-www-form-urlencoded",
      "X-CSRFToken": csrftoken,
    },
    body: `folder=${encodeURIComponent(folder)}`,
  })
    .then((response) => response.json())
    .then((data) => {
      console.log(data);
      if (data.image_paths) {
        data.image_paths.forEach((imagePath, index) => {
          setTimeout(() => {
            const link = document.createElement("a");
            link.href = imagePath;
            link.download = imagePath.split("/").pop();
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
          }, index * 5000);
        });
      }
    })
    .catch((error) => console.error("Error:", error));
}
// fetching image
async function loadImage(folder, filename, containerId) {
  const csrftoken = getCookie("csrftoken");

  await fetch("/get_image/", {
    method: "POST",
    headers: {
      "Content-Type": "application/x-www-form-urlencoded",
      "X-CSRFToken": csrftoken,
    },
    body: `folder=${encodeURIComponent(folder)}&filename=${encodeURIComponent(
      filename
    )}`,
  })
    .then((response) => response.json())
    .then((data) => {
      let container = document.getElementById(containerId);
      if (container) {
        if (data.image_path) {
          if (document.getElementById(filename)) {
            document.getElementById(filename).remove();
          }
          let img = document.createElement("img");
          img.src = data.image_path + `?t=${Date.now()}`;
          img.alt = filename;
          img.id = filename;
          img.classList.add("h-auto", "max-w-full", "mx-auto");
          container.appendChild(img);
        }
      } else {
        console.error("Container not found:", containerId);
      }
    })
    .catch((error) => {
      console.error("Error:", error);
    });
}
// marker option
function insert_marker_options(markers) {
  const dropdownList = document.querySelector("#dropdownmarkers ul");
  dropdownList.innerHTML = "";
  markers.forEach((marker, index) => {
    const li = document.createElement("li");
    li.innerHTML = `
        <div class="flex p-2 rounded hover:bg-gray-100">
          <div class="flex items-center h-5">
            <input
              id="marker-${marker}"
              aria-describedby="marker-text-${index}"
              type="checkbox"
              value="${marker}"
              class="marker-checkbox w-4 h-4 text-blue-600 bg-gray-100 border-gray-300 rounded focus:ring-blue-500 focus:ring-2"
            />
          </div>
          <div class="ms-2 text-sm w-full">
            <label for="marker-${index}" class="font-bold text-gray-900">
              <div>${marker}</div>
            </label>
          </div>
        </div>
      `;
    dropdownList.appendChild(li);
  });
  logSelectedCount();
  dropdownList.addEventListener("change", function (event) {
    if (event.target.type === "checkbox") {
      if (event.target.value === "Select All") {
        handleSelectAll(event.target.checked);
      } else {
        updateSelectAllState();
      }
      logSelectedCount();
    }
  });
}

function handleSelectAll(isChecked) {
  const checkboxes = document.querySelectorAll(".marker-checkbox");
  checkboxes.forEach((checkbox) => {
    checkbox.checked = isChecked;
  });
}

function updateSelectAllState() {
  const checkboxes = document.querySelectorAll(
    '.marker-checkbox:not([value="Select All"])'
  );
  const selectAllCheckbox = document.querySelector(
    '.marker-checkbox[value="Select All"]'
  );
  const allChecked = Array.from(checkboxes).every(
    (checkbox) => checkbox.checked
  );
  selectAllCheckbox.checked = allChecked;
}

function logSelectedCount() {
  const selectedCheckboxes = document.querySelectorAll(
    ".marker-checkbox:checked"
  );
  const selectedMarkers = Array.from(selectedCheckboxes)
    .map((checkbox) => checkbox.value)
    .filter((value) => value !== "Select All");

  // const selectedCount = selectedMarkers.length;

  const selectedbox = document.getElementById("selectedbox");
  selectedbox.textContent = "";
  const markersdiv = document.createElement("div");
  markersdiv.textContent = `${selectedMarkers.join(", ")}`;
  markersdiv.classList.add("indent-8");

  selectedbox.appendChild(markersdiv);
}
